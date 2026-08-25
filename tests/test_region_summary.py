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

from ijump.region_summary import report_average, summarize_by_region

# One cluster per locus: these fixtures predate clusters, and mapping each name
# to itself keeps the pinned values meaningful while the argument changes shape
# (isfinder-annotation 07).
ONE_CLUSTER_EACH = {"IS1": "IS1", "IS2": "IS2"}

GFF_ANN_POS = {
    "tiny_contig": {
        "geneA": ["geneA", "tiny_contig", 800, 840],
        "geneB": ["geneB", "tiny_contig", 841, 900],
    }
}


def _stub_average_depth(chrom, start, stop):
    return 2


def test_summarize_by_region_matches_pinned_golden_output_for_single_hit_path():
    junctions = pd.DataFrame(
        {
            "ID": [1, 2, 3],
            "IS name": ["IS1", "IS2", "IS1"],
            "IS pos": [0, 0, 0],
            "IS_chrom": ["tiny_contig", "tiny_contig", "tiny_contig"],
            "Read name": ["r1", "r2", "r3"],
            "Chrom": ["tiny_contig", "tiny_contig", "tiny_contig"],
            "Position": [820, 850, 780],
            "Orientation": ["l", "l", "l"],
            "Note": ["x", "x", "IS element"],
            "Locus tag": ["a", "a", "a"],
            "Gene": ["g", "g", "g"],
        }
    )

    result = summarize_by_region(junctions, ONE_CLUSTER_EACH, GFF_ANN_POS)

    expected = pd.DataFrame(
        {
            "ann": ["geneA", "geneB"],
            "chrom": ["tiny_contig", "tiny_contig"],
            "start": [800, 841],
            "stop": [840, 900],
            "IS1": [1, 0],
            "IS2": [0, 1],
        },
        index=pd.Index(["geneA", "geneB"], name="ann_id"),
    ).astype(object)

    pdt.assert_frame_equal(result, expected)


def test_report_average_matches_pinned_golden_output_for_single_hit_path():
    sum_by_region = pd.DataFrame(
        {
            "ann": ["geneA", "geneB"],
            "chrom": ["tiny_contig", "tiny_contig"],
            "start": [800, 841],
            "stop": [840, 900],
            "IS1": [1, 0],
            "IS2": [0, 1],
        },
        index=pd.Index(["geneA", "geneB"], name="ann_id"),
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
            "IS Name": ["IS1", "IS2"],
            "Annotation": ["geneA", "geneB"],
            "Chromosome": ["tiny_contig", "tiny_contig"],
            "Start": [800, 841],
            "Stop": [840, 900],
            "Frequency": [4.0, 4.0],
            "Depth": [2, 2],
        },
        index=[0, 3],
    ).astype({"Start": object, "Stop": object})

    pdt.assert_frame_equal(result, expected)


def test_report_average_frequency_is_nan_when_depth_is_zero():
    """Ticket 01: a region can have junction-supporting reads (count > 0)
    but zero total coverage in [start, stop) -- e.g. a fully-deleted
    annotated region. Frequency must come out as NaN there, not inf (which
    plain float division by 0.0 would otherwise produce)."""
    sum_by_region = pd.DataFrame(
        {
            "ann": ["geneA"],
            "chrom": ["tiny_contig"],
            "start": [800],
            "stop": [840],
            "IS1": [1],
            "IS2": [0],
        },
        index=pd.Index(["geneA"], name="ann_id"),
    ).astype(object)

    result = report_average(
        sum_by_region,
        match_lengths=[140, 145, 150],
        read_lengths=3000,
        n_reads_analyzed=20,
        blast_min=10,
        average_depth=lambda chrom, start, stop: 0.0,
    )

    assert result["Depth"].iloc[0] == 0.0
    assert pd.isna(result["Frequency"].iloc[0])


def test_summarize_by_region_accumulates_counts_across_multiple_hits_in_one_region():
    """Spec-based, not characterization -- see module docstring.

    Two junctions for the same IS element (IS1) land in the same annotated
    region (geneA), exercising the "ann_id already registered" branch at
    the count-increment line -- the path this repo's own
    ``pandas-dataframe-append`` rule flagged as untested before this move.
    """
    junctions = pd.DataFrame(
        {
            "ID": [1, 2],
            "IS name": ["IS1", "IS1"],
            "IS pos": [0, 0],
            "IS_chrom": ["tiny_contig", "tiny_contig"],
            "Read name": ["r1", "r2"],
            "Chrom": ["tiny_contig", "tiny_contig"],
            "Position": [810, 830],
            "Orientation": ["l", "l"],
            "Note": ["x", "x"],
            "Locus tag": ["a", "a"],
            "Gene": ["g", "g"],
        }
    )

    result = summarize_by_region(junctions, ONE_CLUSTER_EACH, GFF_ANN_POS)

    expected = pd.DataFrame(
        {
            "ann": ["geneA"],
            "chrom": ["tiny_contig"],
            "start": [800],
            "stop": [840],
            "IS1": [2],
            "IS2": [0],
        },
        index=pd.Index(["geneA"], name="ann_id"),
    ).astype(object)

    pdt.assert_frame_equal(result, expected)


# Two loci of one element plus a second element -- the reference genome's shape,
# where IS17_1, IS17_2 and ISAba12_1 are one ISAba12-like copy and two of its own
# fragments (isfinder-annotation 07).
CLUSTERS = {"IS17_1": "ISAba12", "ISAba12_1": "ISAba12", "ISAlw13_1": "ISAlw13"}


def test_summarize_by_region_emits_one_column_per_cluster():
    """Columns are the elements that jumped, not the loci they were called at.

    A read clipped at one locus of an element cannot be told from one clipped at
    another, so counting them in separate columns splits one insertion's evidence
    across three.
    """
    junctions = _junctions(
        [("IS17_1", 820), ("ISAba12_1", 830), ("ISAlw13_1", 850)],
    )

    result = summarize_by_region(junctions, CLUSTERS, GFF_ANN_POS)

    assert list(result.columns) == ["ann", "chrom", "start", "stop", "ISAba12", "ISAlw13"]


def test_junctions_at_different_loci_of_one_element_accumulate_in_one_column():
    junctions = _junctions([("IS17_1", 820), ("ISAba12_1", 825), ("IS17_1", 830)])

    result = summarize_by_region(junctions, CLUSTERS, GFF_ANN_POS)

    assert result.loc["geneA", "ISAba12"] == 3
    assert result.loc["geneA", "ISAlw13"] == 0


def test_report_average_names_the_cluster_not_the_locus():
    junctions = _junctions([("IS17_1", 820), ("ISAba12_1", 825)])

    summary = summarize_by_region(junctions, CLUSTERS, GFF_ANN_POS)
    report = report_average(summary, [140], 3000, 20, 20, _stub_average_depth)

    assert report["IS Name"].tolist() == ["ISAba12"]
    # One row carrying both junctions, not two rows of one each.
    assert len(report) == 1


def _junctions(is_name_and_position):
    """A junction table carrying just what summarize_by_region reads."""
    return pd.DataFrame(
        {
            "IS name": [is_name for is_name, _ in is_name_and_position],
            "Chrom": ["tiny_contig"] * len(is_name_and_position),
            "Position": [position for _, position in is_name_and_position],
            "Note": ["x"] * len(is_name_and_position),
        }
    )
