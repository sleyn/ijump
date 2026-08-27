"""Origin-spanning elements: the two halves of one IS copy broken by a contig
boundary (isfinder-annotation 05).

The reference genome's case is ``IS17_1`` at the last 144 bases of ``NODE_2``
and ``IS17_2`` at bases 2-77 of the same contig: one ISAba12-like copy the
assembler cut in half. Clustering already unites them; these tests pin the
detection that says *why* they are two rows.
"""

import pandas as pd
import pytest

from ijump import origin_spanning

# NODE_2 of Test/A_baumannii_assembly.fna, and its two IS17 fragments.
NODE_2 = "NODE_2_length_148137_cov_371.33_ID_22131"
NODE_2_LENGTH = 148137


def span(name, start, stop, cluster="ISAba12", contig=NODE_2):
    return origin_spanning.Span(name=name, contig=contig, start=start, stop=stop, cluster=cluster)


def test_a_cluster_touching_both_ends_of_one_contig_spans_the_origin():
    spans = [span("IS17_1", 147994, NODE_2_LENGTH), span("IS17_2", 2, 77)]

    assigned = origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH})

    assert set(assigned) == {"IS17_1", "IS17_2"}
    assert assigned["IS17_1"] == assigned["IS17_2"]


def test_the_element_id_names_the_cluster_it_belongs_to():
    spans = [span("IS17_1", 147994, NODE_2_LENGTH), span("IS17_2", 2, 77)]

    assigned = origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH})

    assert assigned["IS17_1"].startswith("ISAba12")


def test_a_locus_at_a_contig_end_with_no_counterpart_is_not_flagged():
    """The end of a contig is a commonplace; only the pairing makes it a break."""
    spans = [span("IS17_1", 147994, NODE_2_LENGTH), span("ISAba1_1", 112397, 113576)]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_a_locus_at_the_origin_with_no_counterpart_is_not_flagged():
    spans = [span("IS17_2", 2, 77), span("ISAba1_1", 112397, 113576)]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_the_two_halves_have_to_share_a_cluster():
    """Two unrelated elements happening to sit at the two ends of a contig are
    two elements, not one broken in half."""
    spans = [
        span("IS17_1", 147994, NODE_2_LENGTH, cluster="ISAba12"),
        span("ISAlw13_1", 2, 150, cluster="ISAlw13"),
    ]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_the_two_halves_have_to_share_a_contig():
    """A cluster with a copy at the end of one contig and another at the start
    of a second contig is two copies, not one broken in half."""
    other = "NODE_3_length_13486_cov_4113.34_ID_22133"
    spans = [
        span("IS17_1", 147994, NODE_2_LENGTH),
        span("IS17_2", 2, 77, contig=other),
    ]

    assigned = origin_spanning.origin_spanning_elements(
        spans, {NODE_2: NODE_2_LENGTH, other: 13486}
    )

    assert assigned == {}


def test_a_single_locus_covering_a_whole_contig_is_not_flagged():
    """It touches both boundaries, but there is no second span to join to: the
    element is not broken, the contig is just short."""
    short = "NODE_9"
    spans = [span("ISAba12_9", 1, 1039, contig=short)]

    assert origin_spanning.origin_spanning_elements(spans, {short: 1039}) == {}


def test_loci_that_share_a_cluster_but_sit_inside_the_contig_are_not_flagged():
    spans = [span("ISAba11_1", 500, 1600), span("ISAba11_2", 90000, 91100)]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_loci_with_no_cluster_are_left_alone():
    """An empty cluster column says nothing about which rows are one element, so
    it cannot be read as saying they all are."""
    spans = [
        span("IS17_1", 147994, NODE_2_LENGTH, cluster=""),
        span("IS17_2", 2, 77, cluster=""),
    ]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_two_broken_copies_of_one_element_get_separate_ids():
    """Two contigs each broken inside a copy of the same element -- the flag
    groups the halves that belong together, not everything in the cluster."""
    other = "NODE_4"
    spans = [
        span("IS17_1", 147994, NODE_2_LENGTH),
        span("IS17_2", 2, 77),
        span("IS17_3", 4900, 5000, contig=other),
        span("IS17_4", 1, 60, contig=other),
    ]

    assigned = origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH, other: 5000})

    assert assigned["IS17_1"] == assigned["IS17_2"]
    assert assigned["IS17_3"] == assigned["IS17_4"]
    assert assigned["IS17_1"] != assigned["IS17_3"]
    assert len(set(assigned.values())) == 2


def test_a_span_short_of_the_boundary_by_more_than_the_margin_is_not_flagged():
    """The margin exists because alignment ends fray -- ``IS17_2`` starts at base
    2, not 1 -- not to reach loci sitting well inside the contig."""
    spans = [
        span("IS17_1", 147994, NODE_2_LENGTH),
        span("IS17_2", origin_spanning.BOUNDARY_MARGIN + 2, 500),
    ]

    assert origin_spanning.origin_spanning_elements(spans, {NODE_2: NODE_2_LENGTH}) == {}


def test_a_contig_the_lengths_do_not_carry_is_an_error():
    spans = [span("IS17_1", 147994, NODE_2_LENGTH)]

    with pytest.raises(ValueError, match=NODE_2):
        origin_spanning.origin_spanning_elements(spans, {})


# --- the table columns ------------------------------------------------------


def table(rows):
    return pd.DataFrame(rows, columns=["is_name", "contig", "start", "stop", "cluster"])


def test_columns_flag_the_two_halves_and_leave_the_rest_negative():
    flags = origin_spanning.origin_columns(
        table(
            [
                ("IS17_1", NODE_2, "147994", str(NODE_2_LENGTH), "ISAba12"),
                ("IS17_2", NODE_2, "2", "77", "ISAba12"),
                ("ISAba1_1", NODE_2, "112397", "113576", "ISAba1"),
            ]
        ),
        {NODE_2: NODE_2_LENGTH},
    )

    assert list(flags["wraps_origin"]) == ["yes", "yes", "no"]
    assert flags["element_id"][0] == flags["element_id"][1]
    assert flags["element_id"][2] == ""


def test_columns_keep_both_coordinate_rows_untouched():
    """The two spans are never joined into one ``start > stop`` row: the
    boundary search and the per-region overlap logic both assume a row's
    start precedes its stop."""
    rows = table(
        [
            ("IS17_1", NODE_2, "147994", str(NODE_2_LENGTH), "ISAba12"),
            ("IS17_2", NODE_2, "2", "77", "ISAba12"),
        ]
    )

    flags = origin_spanning.origin_columns(rows, {NODE_2: NODE_2_LENGTH})

    assert len(flags) == len(rows)
    assert list(rows["start"]) == ["147994", "2"]
    assert list(rows["stop"]) == [str(NODE_2_LENGTH), "77"]
