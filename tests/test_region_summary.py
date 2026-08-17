"""Tests for ticket 13: ``region_summary.summarize_by_region``/``report_average``.

The single-hit tests are characterization tests: the fixtures were run against
today's (pre-move) ``ISClipped.summary_junctions_by_region``/``report_average``
methods (the ``.append()``-based implementation) to pin the expected output
before the move -- they document today's behaviour, not whether it's correct.
Every region here receives at most one junction total, matching the ticket's
note that ``isclipped.py:660-662``'s "already registered" branch is never hit
in practice.

The multi-hit test is spec-based, not characterization: no historical baseline
exists for two-or-more junctions landing in the same annotated region (nothing
exercises that path today), so it asserts against the new
``pd.concat``-based implementation directly.
"""
import pandas as pd
import pandas.testing as pdt

from ijump.region_summary import summarize_by_region, report_average

IS_COORDS = {'IS1': ['tiny_contig', '1', '10'], 'IS2': ['tiny_contig', '1', '10']}

GFF_ANN_POS = {
    'tiny_contig': {
        'geneA': ['geneA', 'tiny_contig', 800, 840],
        'geneB': ['geneB', 'tiny_contig', 841, 900],
    }
}


def _stub_average_depth(chrom, start, stop):
    return 2


def test_summarize_by_region_matches_pinned_golden_output_for_single_hit_path():
    junctions = pd.DataFrame({
        'ID': [1, 2, 3],
        'IS name': ['IS1', 'IS2', 'IS1'],
        'IS pos': [0, 0, 0],
        'IS_chrom': ['tiny_contig', 'tiny_contig', 'tiny_contig'],
        'Read name': ['r1', 'r2', 'r3'],
        'Chrom': ['tiny_contig', 'tiny_contig', 'tiny_contig'],
        'Position': [820, 850, 780],
        'Orientation': ['l', 'l', 'l'],
        'Note': ['x', 'x', 'IS element'],
        'Locus tag': ['a', 'a', 'a'],
        'Gene': ['g', 'g', 'g'],
    })

    result = summarize_by_region(junctions, IS_COORDS, GFF_ANN_POS)

    expected = pd.DataFrame(
        {
            'ann': ['geneA', 'geneB'],
            'chrom': ['tiny_contig', 'tiny_contig'],
            'start': [800, 841],
            'stop': [840, 900],
            'IS1': [1, 0],
            'IS2': [0, 1],
        },
        index=pd.Index(['geneA', 'geneB'], name='ann_id'),
    ).astype(object)

    pdt.assert_frame_equal(result, expected)


def test_report_average_matches_pinned_golden_output_for_single_hit_path():
    sum_by_region = pd.DataFrame(
        {
            'ann': ['geneA', 'geneB'],
            'chrom': ['tiny_contig', 'tiny_contig'],
            'start': [800, 841],
            'stop': [840, 900],
            'IS1': [1, 0],
            'IS2': [0, 1],
        },
        index=pd.Index(['geneA', 'geneB'], name='ann_id'),
    ).astype(object)

    result = report_average(
        sum_by_region,
        match_lengths=[140, 145, 150],
        read_lengths=3000,
        n_reads_analyzed=20,
        blast_min=10,
        average_depth=_stub_average_depth,
    )

    expected = pd.DataFrame(
        {
            'IS Name': ['IS1', 'IS2'],
            'Annotation': ['geneA', 'geneB'],
            'Chromosome': ['tiny_contig', 'tiny_contig'],
            'Start': [800, 841],
            'Stop': [840, 900],
            'Frequency': [4.0, 4.0],
            'Depth': [2, 2],
        },
        index=[0, 3],
    ).astype({'Start': object, 'Stop': object})

    pdt.assert_frame_equal(result, expected)


def test_summarize_by_region_accumulates_counts_across_multiple_hits_in_one_region():
    """Spec-based, not characterization -- see module docstring.

    Two junctions for the same IS element (IS1) land in the same annotated
    region (geneA), exercising the "ann_id already registered" branch at
    the count-increment line -- the path this repo's own
    ``pandas-dataframe-append`` rule flagged as untested before this move.
    """
    junctions = pd.DataFrame({
        'ID': [1, 2],
        'IS name': ['IS1', 'IS1'],
        'IS pos': [0, 0],
        'IS_chrom': ['tiny_contig', 'tiny_contig'],
        'Read name': ['r1', 'r2'],
        'Chrom': ['tiny_contig', 'tiny_contig'],
        'Position': [810, 830],
        'Orientation': ['l', 'l'],
        'Note': ['x', 'x'],
        'Locus tag': ['a', 'a'],
        'Gene': ['g', 'g'],
    })

    result = summarize_by_region(junctions, IS_COORDS, GFF_ANN_POS)

    expected = pd.DataFrame(
        {
            'ann': ['geneA'],
            'chrom': ['tiny_contig'],
            'start': [800],
            'stop': [840],
            'IS1': [2],
            'IS2': [0],
        },
        index=pd.Index(['geneA'], name='ann_id'),
    ).astype(object)

    pdt.assert_frame_equal(result, expected)
