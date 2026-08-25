"""Characterization test for ticket 06: ``junction_pairing.find_pairs``.

Test/test_find_pair.py (gitignored, not run by anything) pasted its own copy
of _find_pair rather than importing the shipped one, and had already drifted
from it -- the copy builds 'Count_l'/'Count_r' columns while the real
function in junction_pairing.py uses 'Count_mapped_to_IS_l'/'_r'. It also
had no assertions: it called the function, bound the result to a variable,
and printed 'OK' unconditionally.

This test imports the real function and pins its output on the same
fixture (captured from actual data) as a golden value. It is characterization,
not specification -- it documents today's behaviour, not whether that
behaviour is correct.
"""

import numpy as np
import pandas as pd
import pandas.testing as pdt

from ijump.junction_pairing import find_pairs

POS_L = np.array([311753, 311755, 311759, 311773, 311788, 992367, 3352696, 3790412])
POS_L_COUNT = np.array([1, 190, 1, 1, 1, 1, 1, 1])
POS_R = np.array([311758, 311760, 311761, 311762, 311763, 318981, 586311, 1190170, 3790421])
POS_R_COUNT = np.array([1, 1, 2, 1, 154, 2, 1, 1, 2])
CHROM_LEN = 3909467
MAX_IS_DUP_LEN = 20
CHROM = "NODE_1_length_3909467_cov_533.478_ID_22129"


def test_find_pair_matches_pinned_golden_output():
    pairs_df = find_pairs(
        POS_L.copy(),
        POS_R.copy(),
        POS_L_COUNT.copy(),
        POS_R_COUNT.copy(),
        CHROM_LEN,
        MAX_IS_DUP_LEN,
        CHROM,
    )

    expected = pd.DataFrame(
        {
            "Position_l": [
                311755,
                311773,
                311759,
                311753,
                311788,
                992367,
                3352696,
                3790412,
                0,
                0,
                0,
                0,
            ],
            "Position_r": [
                311763,
                311758,
                311760,
                311762,
                0,
                0,
                0,
                3790421,
                311761,
                318981,
                586311,
                1190170,
            ],
            "Count_mapped_to_IS_l": [190, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
            "Count_mapped_to_IS_r": [154, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1],
            "Chrom": [CHROM] * 12,
        }
    )

    pdt.assert_frame_equal(pairs_df, expected)


def test_find_pair_accepts_read_only_input_arrays():
    """pandas 3 hands out read-only arrays, and precise mode feeds them straight in.

    ``isclipped.search_insert_pos`` builds every argument here with
    ``Series.to_numpy()``, which under pandas 3 returns a non-writeable view.
    ``find_pairs`` sorts its clusters by writing back into ``pos_l``,
    ``pos_l_count`` and ``closeness_matrix``, so a non-writeable input aborted
    the whole precise run with "assignment destination is read-only" -- before
    a single output file was written.
    """
    read_only_args = []
    for array in (POS_L, POS_R, POS_L_COUNT, POS_R_COUNT):
        copy = array.copy()
        copy.flags.writeable = False
        read_only_args.append(copy)

    pairs_df = find_pairs(*read_only_args, CHROM_LEN, MAX_IS_DUP_LEN, CHROM)

    expected = find_pairs(
        POS_L.copy(),
        POS_R.copy(),
        POS_L_COUNT.copy(),
        POS_R_COUNT.copy(),
        CHROM_LEN,
        MAX_IS_DUP_LEN,
        CHROM,
    )
    pdt.assert_frame_equal(pairs_df, expected)


def test_find_pair_leaves_its_input_arrays_untouched():
    """The caller's arrays are not scratch space -- the cluster sort must not reorder them."""
    pos_l, pos_r = POS_L.copy(), POS_R.copy()
    pos_l_count, pos_r_count = POS_L_COUNT.copy(), POS_R_COUNT.copy()

    find_pairs(pos_l, pos_r, pos_l_count, pos_r_count, CHROM_LEN, MAX_IS_DUP_LEN, CHROM)

    np.testing.assert_array_equal(pos_l, POS_L)
    np.testing.assert_array_equal(pos_r, POS_R)
    np.testing.assert_array_equal(pos_l_count, POS_L_COUNT)
    np.testing.assert_array_equal(pos_r_count, POS_R_COUNT)


# One-sided input: every junction on the contig points the same way, so there is
# nothing to pair and ``find_pairs`` returns early. That exit has one arm per
# side (junction-pairing-orphans 01).
ONE_SIDED_POS = np.array([451, 145216, 145218])
ONE_SIDED_COUNT = np.array([1, 7, 3])
EMPTY = np.array([], dtype=np.int64)


def test_right_only_input_is_written_as_right_orphans():
    """A right junction belongs in the right columns.

    The right-only arm used to write position and count into ``Position_l`` and
    ``Count_mapped_to_IS_l``, mirroring the left-only arm rather than the main
    path's own right-orphan writer. Frequency estimation reads a right
    junction's evidence out of the right columns, found zeros there, and
    reported 0.0 -- so the slip silently zeroed real detections.
    """
    pairs_df = find_pairs(
        EMPTY, ONE_SIDED_POS, EMPTY, ONE_SIDED_COUNT, CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )

    expected = pd.DataFrame(
        {
            "Position_l": [0, 0, 0],
            "Position_r": ONE_SIDED_POS,
            "Count_mapped_to_IS_l": [0, 0, 0],
            "Count_mapped_to_IS_r": ONE_SIDED_COUNT,
            "Chrom": CHROM,
        }
    )
    pdt.assert_frame_equal(pairs_df, expected, check_dtype=False)


def test_left_only_input_is_written_as_left_orphans():
    """The left-only arm was already right; pinned so the fix to its sibling
    cannot be applied to both."""
    pairs_df = find_pairs(
        ONE_SIDED_POS, EMPTY, ONE_SIDED_COUNT, EMPTY, CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )

    expected = pd.DataFrame(
        {
            "Position_l": ONE_SIDED_POS,
            "Position_r": [0, 0, 0],
            "Count_mapped_to_IS_l": ONE_SIDED_COUNT,
            "Count_mapped_to_IS_r": [0, 0, 0],
            "Chrom": CHROM,
        }
    )
    pdt.assert_frame_equal(pairs_df, expected, check_dtype=False)


def test_one_sided_exit_writes_an_orphan_the_way_the_main_path_does():
    """The two exits have to agree: an unpaired right junction looks the same
    whether or not the contig happened to carry left junctions elsewhere.
    Both build their rows through ``_right_orphan_row``, so this holds that
    shared shape rather than two copies of it.

    Driven through the main path by adding a left junction far enough away that
    nothing pairs, then comparing the right rows the two exits produce.
    """
    far_left = np.array([2_000_000])
    main_path = find_pairs(
        far_left,
        ONE_SIDED_POS,
        np.array([1]),
        ONE_SIDED_COUNT,
        CHROM_LEN,
        MAX_IS_DUP_LEN,
        CHROM,
    )
    early_return = find_pairs(
        EMPTY, ONE_SIDED_POS, EMPTY, ONE_SIDED_COUNT, CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )

    right_rows = main_path.query("Position_l == 0").reset_index(drop=True)
    pdt.assert_frame_equal(right_rows, early_return, check_dtype=False)
