"""Characterization test for ticket 09: ``frequency_estimation.estimate_frequencies``.

The fixture below is hand-built (synthetic), not observed real data -- there is no
pre-existing captured-real-data fixture for this function, and the tiny end-to-end
fixture in tests/fixtures/ never reaches this code path (it exits at the
no-clipped-reads path, per ticket 03's note). The expected values were pinned by
running today's (pre-move) ``ISClipped.assess_isel_freq`` body against this same
fixture. This is characterization, not specification -- it documents today's
behaviour, not whether that behaviour is correct.
"""
import pandas as pd
import pandas.testing as pdt

from frequency_estimation import estimate_frequencies

PAIRS_DF = pd.DataFrame({
    'Position_l': [100, 150, 0],
    'Position_r': [200, 0, 250],
    'Count_mapped_to_IS_l': [5, 3, 0],
    'Count_mapped_to_IS_r': [4, 0, 2],
    'Chrom': ['chrA', 'chrA', 'chrA'],
    'IS_name': ['IS1', 'IS1', 'IS1'],
})

CLIPPED_READS_BWRD = pd.DataFrame({
    'clip_position': ['left', 'left', 'left', 'right', 'right', 'right'],
    'junction_in_read': [100, 100, 150, 200, 200, 250],
    'IS_chrom': ['chrA', 'chrA', 'chrA', 'chrA', 'chrA', 'chrA'],
    'Read name': ['r1', 'r2', 'r3', 'r4', 'r5', 'r6'],
})

UNCLIPPED_DEPTH = {'chrA': {100: 10, 150: 6, 200: 8, 250: 3}}
CL_READ_COV_OVERLAP = {'chrA': {100: 2, 150: 1, 200: 3, 250: 0}}
MATCH_LENGTHS = [140, 145, 150]
READ_LENGTHS = 3000
N_READS_ANALYZED = 20


def test_estimate_frequencies_matches_pinned_golden_output():
    result = estimate_frequencies(
        PAIRS_DF.copy(), CLIPPED_READS_BWRD, UNCLIPPED_DEPTH, CL_READ_COV_OVERLAP,
        MATCH_LENGTHS, READ_LENGTHS, N_READS_ANALYZED,
    )

    expected = pd.DataFrame({
        'Position_l': [100, 150, 0],
        'Position_r': [200, 0, 250],
        'Count_mapped_to_IS_l': [5, 3, 0],
        'Count_mapped_to_IS_r': [4, 0, 2],
        'Chrom': ['chrA', 'chrA', 'chrA'],
        'IS_name': ['IS1', 'IS1', 'IS1'],
        'N_unclipped_l': [10, 6, 0],
        'N_clipped_l': [2.0, 1.0, 0.0],
        'N_unclipped_r': [8, 0, 3],
        'N_clipped_r': [2.0, 0.0, 1.0],
        'N_overlap_l': [2, 1, 0],
        'N_overlap_r': [3, 0, 0],
        'N_clipped_l_correction': [28.000000000000007, 14.000000000000004, 0.0],
        'N_clipped_r_correction': [28.000000000000007, 0.0, 14.000000000000004],
        'N_overlap_l_correction': [28.000000000000007, 14.000000000000004, 0.0],
        'N_overlap_r_correction': [42.00000000000001, 0.0, 0.0],
        'N_clipped_l_corrected': [30.000000000000007, 15.000000000000004, 0.0],
        'N_overlap_l_corrected': [30.000000000000007, 15.000000000000004, 0.0],
        'N_clipped_r_corrected': [30.000000000000007, 0.0, 15.000000000000004],
        'N_overlap_r_corrected': [45.00000000000001, 0.0, 0.0],
        'N_overlap_formula_l': [0.0, 15.000000000000004, 0.0],
        'N_overlap_formula_r': [15.0, 0.0, 0.0],
        'Frequency_l': [0.7481296758104738, 0.4155124653739612, 0.0],
        'Frequency_r': [0.5649717514124294, 0.0, 0.8287292817679558],
        'Frequency': [0.6565507136114517, 0.4155124653739612, 0.8287292817679558],
        'Depth': [46.50000000000001, 36.00000000000001, 18.000000000000004],
    })

    pdt.assert_frame_equal(result, expected)


def test_estimate_frequencies_does_not_mutate_input_pairs_df():
    original = PAIRS_DF.copy()

    estimate_frequencies(
        PAIRS_DF, CLIPPED_READS_BWRD, UNCLIPPED_DEPTH, CL_READ_COV_OVERLAP,
        MATCH_LENGTHS, READ_LENGTHS, N_READS_ANALYZED,
    )

    pdt.assert_frame_equal(PAIRS_DF, original)
