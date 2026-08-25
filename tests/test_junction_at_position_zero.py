"""A junction at the first base of a contig is a position, not an absence.

Positions inside the pipeline are 0-based, so `0` is the first base of a contig.
It was also the marker for "this junction has no partner", and the two meanings
collided: a real junction at position 0 read as missing everywhere the marker was
tested for (junction-pairing-orphans 02).

Absence is now `NO_JUNCTION`, outside the coordinate space. The output file is
unchanged -- it is 1-based, where 0 is not a position and so is an unambiguous
"absent" already.
"""

import numpy as np
import pytest
from fake_clipped_read import FakeRead

from ijump.clipped_read_search import _clboundaries
from ijump.isclipped import convert_zero_one_base, interpos_distance, keep_pair
from ijump.junction_pairing import NO_JUNCTION, find_pairs

CHROM = "contig_1"
CHROM_LEN = 3909467
MAX_IS_DUP_LEN = 20
EMPTY = np.array([], dtype=np.int64)


def test_a_read_clipped_at_the_first_base_of_a_contig_yields_position_zero():
    """Reachability, established rather than assumed.

    ``get_reference_positions`` is 0-based, so a left-clipped read whose aligned
    part starts at the contig's first base reports its junction at 0. This is the
    origin-spanning shape: the assembler's break lands inside an IS copy, leaving
    a fragment at the very start of the contig (isfinder-annotation 05). The
    reference genome comes within one base of it -- ``IS17_2`` starts at base 2
    of ``NODE_2``.
    """
    read = FakeRead("r1", "3S5M", [None, None, None, 0, 1, 2, 3, 4], "AAACCCCC")

    assert _clboundaries(read) == [[1, 3, "left", 0]]


def test_zero_is_a_position_on_the_one_sided_exit():
    """The exit taken when nothing on the contig points the other way."""
    pairs_df = find_pairs(
        np.array([0]), EMPTY, np.array([7]), EMPTY, CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )

    assert len(pairs_df) == 1
    assert pairs_df.iloc[0]["Position_l"] == 0
    assert pairs_df.iloc[0]["Position_r"] == NO_JUNCTION
    assert pairs_df.iloc[0]["Count_mapped_to_IS_l"] == 7


def test_zero_is_a_position_on_the_main_path():
    """The same junction, on the exit taken when the contig does carry junctions
    on the other side. It used to be dropped here and kept there."""
    pairs_df = find_pairs(
        np.array([0]),
        np.array([500_000]),
        np.array([7]),
        np.array([3]),
        CHROM_LEN,
        MAX_IS_DUP_LEN,
        CHROM,
    )

    left = pairs_df.query("Position_l == 0")
    assert len(left) == 1
    assert left.iloc[0]["Position_r"] == NO_JUNCTION
    assert left.iloc[0]["Count_mapped_to_IS_l"] == 7


def test_both_exits_write_the_same_row_for_a_junction_at_zero():
    one_sided = find_pairs(
        np.array([0]), EMPTY, np.array([7]), EMPTY, CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )
    main_path = find_pairs(
        np.array([0]),
        np.array([500_000]),
        np.array([7]),
        np.array([3]),
        CHROM_LEN,
        MAX_IS_DUP_LEN,
        CHROM,
    )

    from_main_path = main_path.query("Position_l == 0").reset_index(drop=True)
    assert from_main_path.to_dict("records") == one_sided.to_dict("records")


def test_the_main_path_still_drops_the_rows_it_never_wrote():
    """The blank rows the pre-sized frame carries are not junctions at 0.

    They used to be cleared by dropping every row with both positions zero,
    which is what swallowed the real thing.
    """
    pos_l = np.array([1000, 2000])
    pos_r = np.array([1005, 3000, 4000])

    pairs_df = find_pairs(
        pos_l, pos_r, np.array([1, 1]), np.array([1, 1, 1]), CHROM_LEN, MAX_IS_DUP_LEN, CHROM
    )

    # Two left junctions and one right orphan each, no blanks: the frame is
    # sized for pos_l.size + pos_r.size = 5.
    assert len(pairs_df) == 4
    assert not (
        (pairs_df["Position_l"] == NO_JUNCTION) & (pairs_df["Position_r"] == NO_JUNCTION)
    ).any()


@pytest.mark.parametrize(
    "coordinates,expected",
    [
        # A real 0-based position becomes its 1-based self.
        ([0, 4], [1, 5]),
        # Absence stays absent: 0 in a 1-based file is not a position.
        ([NO_JUNCTION], [0]),
    ],
)
def test_base_conversion_separates_absence_from_the_first_base(coordinates, expected):
    assert convert_zero_one_base(coordinates) == expected


def test_interpos_distance_reads_zero_as_a_position():
    """A pair spanning the contig's first base has a real span; only an absent
    side falls back to the default."""
    assert interpos_distance(0, 8) == 13
    assert interpos_distance(0, NO_JUNCTION) == 5
    assert interpos_distance(NO_JUNCTION, 8) == 5


def test_keep_pair_does_not_match_a_region_on_an_absent_side():
    """Regions can start at 0, so an absent side encoded as 0 fell inside them
    and kept pairs on evidence that was not there."""
    starts, ends = np.array([0]), np.array([100])

    assert keep_pair([0, NO_JUNCTION], starts, ends) is True
    assert keep_pair([NO_JUNCTION, 500], starts, ends) is False
