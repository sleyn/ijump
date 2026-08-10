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

from junction_pairing import find_pairs

POS_L = np.array([311753, 311755, 311759, 311773, 311788, 992367, 3352696, 3790412])
POS_L_COUNT = np.array([1, 190, 1, 1, 1, 1, 1, 1])
POS_R = np.array([311758, 311760, 311761, 311762, 311763, 318981, 586311, 1190170, 3790421])
POS_R_COUNT = np.array([1, 1, 2, 1, 154, 2, 1, 1, 2])
CHROM_LEN = 3909467
MAX_IS_DUP_LEN = 20
CHROM = 'NODE_1_length_3909467_cov_533.478_ID_22129'


def test_find_pair_matches_pinned_golden_output():
    pairs_df = find_pairs(
        POS_L.copy(), POS_R.copy(), POS_L_COUNT.copy(), POS_R_COUNT.copy(),
        CHROM_LEN, MAX_IS_DUP_LEN, CHROM,
    )

    expected = pd.DataFrame({
        'Position_l': [311755, 311773, 311759, 311753, 311788, 992367, 3352696, 3790412, 0, 0, 0, 0],
        'Position_r': [311763, 311758, 311760, 311762, 0, 0, 0, 3790421, 311761, 318981, 586311, 1190170],
        'Count_mapped_to_IS_l': [190, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'Count_mapped_to_IS_r': [154, 1, 1, 1, 0, 0, 0, 2, 2, 2, 1, 1],
        'Chrom': [CHROM] * 12,
    })

    pdt.assert_frame_equal(pairs_df, expected)
