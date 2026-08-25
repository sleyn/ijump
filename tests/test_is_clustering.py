"""Similarity clustering of the called IS loci (isfinder-annotation 04).

Two loci belong to the same mobile element when they align at >=95% identity
over >=80% of the *shorter* one, and clusters are closed under single linkage --
so a fragment reaches a parent it shares no alignment with, through a sibling
that aligns to both.

The numbers in these tests are the measured all-vs-all ``blastn`` of the 13 loci
called for ``Test/A_baumannii_assembly.fna`` (recorded in the spec), so the cases
the ticket names are exercised here without the 4 MB genome.
"""

import logging

import pandas as pd
import pytest

from ijump import is_clustering

# The 13 loci of Test/A_baumannii_assembly.fna, as (name, contig, start, stop).
REFERENCE_LOCI = [
    ("ISAba18_1", "NODE_2", 93700, 95008),
    ("ISAba1_1", "NODE_2", 112397, 113576),
    ("ISAba11_1", "NODE_1", 1719213, 1720313),
    ("ISAba11_2", "NODE_1", 952634, 953734),
    ("ISAba11_3", "NODE_2", 134383, 135483),
    ("ISAba53_1", "NODE_1", 708843, 709881),
    ("ISAba12_1", "NODE_1", 2983841, 2984879),
    ("ISVsa3_1", "NODE_2", 110032, 111008),
    ("ISAcsp3_1", "NODE_1", 2980551, 2981283),
    ("ISAlw32_1", "NODE_1", 2980018, 2980550),
    ("IS17_1", "NODE_2", 147994, 148137),
    ("ISAlw13_1", "NODE_3", 553, 701),
    ("IS17_2", "NODE_2", 2, 77),
]

# The pairs that align at all. Every other pair has no hit.
REFERENCE_HITS = [
    ("ISAba11_1", "ISAba11_2", 100.0, 1101),
    ("ISAba11_2", "ISAba11_3", 100.0, 1101),
    ("ISAba11_1", "ISAba11_3", 100.0, 1101),
    ("IS17_2", "ISAba12_1", 100.0, 76),
    ("IS17_1", "ISAba12_1", 97.9, 144),
    ("ISAba53_1", "ISAba12_1", 83.0, 1045),
]


def reference_loci():
    return [
        is_clustering.Locus(name=name, contig=contig, start=start, length=stop - start + 1)
        for name, contig, start, stop in REFERENCE_LOCI
    ]


def reference_hits():
    """Every recorded hit, in both directions -- blastn reports each pair twice."""
    hits = []
    for query, subject, identity, length in REFERENCE_HITS:
        hits.append(is_clustering.Hit(query, subject, identity, length))
        hits.append(is_clustering.Hit(subject, query, identity, length))
    return hits


def clusters_by_member(clusters):
    """``{locus name: frozenset of its cluster's members}``, for set comparisons."""
    return {
        locus.name: frozenset(member.name for member in cluster.members)
        for cluster in clusters
        for locus in cluster.members
    }


# --- linkage ---------------------------------------------------------------


def test_is17_fragments_join_isaba12_and_isaba53_stays_out():
    """The case the whole change exists for.

    ``IS17_1`` and ``IS17_2`` are the two opposite ends of one ISAba12-like
    element and do not align to each other at all -- they reach each other only
    through ``ISAba12_1``. ``ISAba53_1`` aligns to ``ISAba12_1`` over its whole
    length but at 83%, which a clipped read tells apart, so it stays separate.
    """
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())
    members = clusters_by_member(clusters)

    assert members["ISAba12_1"] == frozenset({"IS17_1", "IS17_2", "ISAba12_1"})
    assert members["ISAba53_1"] == frozenset({"ISAba53_1"})


def test_identical_isaba11_copies_share_one_cluster():
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())
    members = clusters_by_member(clusters)

    assert members["ISAba11_1"] == frozenset({"ISAba11_1", "ISAba11_2", "ISAba11_3"})


def test_every_locus_lands_in_exactly_one_cluster():
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())

    placed = [locus.name for cluster in clusters for locus in cluster.members]
    assert sorted(placed) == sorted(name for name, _, _, _ in REFERENCE_LOCI)


def test_a_locus_with_no_hit_at_all_is_its_own_cluster():
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())
    members = clusters_by_member(clusters)

    assert members["ISAba18_1"] == frozenset({"ISAba18_1"})


def test_coverage_is_measured_on_the_shorter_locus():
    """76 bp of alignment is 7% of ``ISAba12_1`` but 100% of the remnant.

    Measuring on the longer locus would drop the remnant, which is the wrong
    answer: a read clipped at a 76 bp remnant cannot be told from one clipped at
    the full copy.
    """
    loci = [
        is_clustering.Locus("remnant_1", "NODE_1", 1, 76),
        is_clustering.Locus("parent_1", "NODE_1", 500, 1039),
    ]
    hits = [
        is_clustering.Hit("remnant_1", "parent_1", 100.0, 76),
        is_clustering.Hit("parent_1", "remnant_1", 100.0, 76),
    ]

    clusters = is_clustering.cluster_loci(loci, hits)

    assert len(clusters) == 1


@pytest.mark.parametrize(
    "identity,aligned,linked",
    [
        (100.0, 1039, True),
        # Just under either threshold and the pair is not linked.
        (94.9, 1039, False),
        (100.0, 832, True),  # 80% of 1039 is 831.2
        (100.0, 831, False),
    ],
)
def test_both_thresholds_gate_the_link(identity, aligned, linked):
    loci = [
        is_clustering.Locus("a_1", "NODE_1", 1, 1039),
        is_clustering.Locus("b_1", "NODE_1", 5000, 1039),
    ]
    hits = [is_clustering.Hit("a_1", "b_1", identity, aligned)]

    clusters = is_clustering.cluster_loci(loci, hits)

    assert (len(clusters) == 1) is linked


def test_thresholds_are_configurable():
    """``ISAba53_1`` and ``ISAba12_1`` merge once the identity floor drops
    below their 83%, which is what the flags are for."""
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits(), identity=80.0)
    members = clusters_by_member(clusters)

    assert "ISAba53_1" in members["ISAba12_1"]


def test_a_hit_of_a_locus_to_itself_does_not_link_anything():
    loci = [is_clustering.Locus("a_1", "NODE_1", 1, 1039)]
    hits = [is_clustering.Hit("a_1", "a_1", 100.0, 1039)]

    clusters = is_clustering.cluster_loci(loci, hits)

    assert len(clusters) == 1
    assert [locus.name for locus in clusters[0].members] == ["a_1"]


# --- chaining ---------------------------------------------------------------


def test_chained_pairs_are_reported_by_name():
    """Single linkage lets a pair that fails the threshold share a cluster.

    That is wanted for ``IS17_1``/``IS17_2`` and unwanted when a conserved
    stretch bridges two distinct elements -- and nothing in the alignment tells
    the two apart, so every such pair is surfaced for the operator to judge.
    """
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())
    (isaba12,) = [c for c in clusters if any(m.name == "ISAba12_1" for m in c.members)]

    assert isaba12.chained_pairs == [("IS17_1", "IS17_2")]


def test_a_cluster_whose_every_pair_passes_reports_no_chaining():
    clusters = is_clustering.cluster_loci(reference_loci(), reference_hits())
    (isaba11,) = [c for c in clusters if any(m.name == "ISAba11_1" for m in c.members)]

    assert isaba11.chained_pairs == []


def test_chained_pairs_are_logged_with_both_locus_names(caplog):
    with caplog.at_level(logging.WARNING):
        is_clustering.cluster_column(reference_loci(), reference_hits())

    warnings = [record.getMessage() for record in caplog.records]
    assert any("IS17_1" in message and "IS17_2" in message for message in warnings)


# --- naming -----------------------------------------------------------------


@pytest.mark.parametrize(
    "is_name,expected",
    [
        ("ISAba12_1", "ISAba12"),
        ("ISAba11_12", "ISAba11"),
        # 11 database names carry an underscore of their own; only the copy
        # suffix the parser appended comes off.
        ("ISBj2_B_1", "ISBj2_B"),
        # Nothing that is not a copy suffix is stripped.
        ("ISAba12", "ISAba12"),
        ("IS1247a_ISPpu12", "IS1247a_ISPpu12"),
    ],
)
def test_base_is_name_strips_only_the_copy_suffix(is_name, expected):
    assert is_clustering.base_is_name(is_name) == expected


def test_cluster_takes_its_longest_members_name():
    """``{IS17_1, IS17_2, ISAba12_1}`` is one ISAba12-like element with two of
    its own fragments, so the full copy names it -- not the fragment."""
    column = is_clustering.cluster_column(reference_loci(), reference_hits())

    assert column["ISAba12_1"] == "ISAba12"
    assert column["IS17_1"] == "ISAba12"
    assert column["IS17_2"] == "ISAba12"


def test_identical_copies_are_named_by_the_first_coordinate():
    """All three ``ISAba11`` copies are 1101 bp, so length settles nothing and
    the tie falls to the coordinate -- a rule, so the name does not depend on
    the order rows happen to arrive in."""
    column = is_clustering.cluster_column(reference_loci(), reference_hits())

    assert column["ISAba11_1"] == "ISAba11"
    assert column["ISAba11_3"] == "ISAba11"


def test_a_lone_cluster_keeps_the_bare_name():
    column = is_clustering.cluster_column(reference_loci(), reference_hits())

    assert column["ISAba18_1"] == "ISAba18"
    assert column["ISAba53_1"] == "ISAba53"


def test_colliding_clusters_are_suffixed_by_descending_size():
    """The database is sparse relative to real IS diversity, so two clusters
    ~85% apart can share a nearest database neighbour and both claim its name.
    Suffixes appear only then, largest cluster first."""
    loci = [
        is_clustering.Locus("ISAba12_1", "NODE_1", 1000, 1039),
        is_clustering.Locus("ISAba12_2", "NODE_1", 5000, 1039),
        is_clustering.Locus("ISAba12_3", "NODE_1", 9000, 1039),
    ]
    # _1 and _2 are one element; _3 is a different one that happened to hit the
    # same database entry.
    hits = [
        is_clustering.Hit("ISAba12_1", "ISAba12_2", 100.0, 1039),
        is_clustering.Hit("ISAba12_2", "ISAba12_1", 100.0, 1039),
    ]

    column = is_clustering.cluster_column(loci, hits)

    assert column["ISAba12_1"] == "ISAba12.a"
    assert column["ISAba12_2"] == "ISAba12.a"
    assert column["ISAba12_3"] == "ISAba12.b"


def test_equal_sized_colliding_clusters_are_suffixed_by_coordinate():
    loci = [
        is_clustering.Locus("ISAba12_1", "NODE_1", 9000, 1039),
        is_clustering.Locus("ISAba12_2", "NODE_1", 1000, 1039),
    ]

    column = is_clustering.cluster_column(loci, [])

    assert column["ISAba12_2"] == "ISAba12.a"
    assert column["ISAba12_1"] == "ISAba12.b"


def test_cluster_names_are_unique():
    """A duplicate is a correctness bug, not a cosmetic one: ``is_coords`` is a
    dict keyed on this name, so one locus would silently overwrite another."""
    # Every locus a separate cluster, all claiming the same database name.
    loci = [is_clustering.Locus(f"ISAba12_{n}", "NODE_1", n * 5000, 1039) for n in range(1, 31)]

    column = is_clustering.cluster_column(loci, [])

    assert len(set(column.values())) == len(loci)


def test_more_collisions_than_letters_still_yield_unique_names():
    """26 colliding clusters exhaust a single letter; the 27th must not wrap
    back onto the first."""
    loci = [is_clustering.Locus(f"ISAba12_{n}", "NODE_1", n * 5000, 1039) for n in range(1, 31)]

    column = is_clustering.cluster_column(loci, [])

    assert column["ISAba12_1"] == "ISAba12.a"
    assert column["ISAba12_26"] == "ISAba12.z"
    assert column["ISAba12_27"] == "ISAba12.aa"


# --- table input ------------------------------------------------------------


def test_a_table_with_two_rows_of_one_name_is_rejected():
    """Everything downstream keys on ``is_name`` -- linkage here, ``is_coords``
    later -- so a duplicate would drop a row instead of failing."""
    table = pd.DataFrame(
        [
            ["ISAba12_1", "NODE_1", "100", "1138"],
            ["ISAba12_1", "NODE_1", "5000", "6038"],
        ],
        columns=["is_name", "contig", "start", "stop"],
    )

    with pytest.raises(ValueError, match="ISAba12_1"):
        is_clustering.loci_from_table(table)


def test_table_rows_become_loci_with_inclusive_lengths():
    table = pd.DataFrame(
        [["ISAba12_1", "NODE_1", "100", "1138"]],
        columns=["is_name", "contig", "start", "stop"],
    )

    (locus,) = is_clustering.loci_from_table(table)

    assert locus == is_clustering.Locus("ISAba12_1", "NODE_1", 100, 1039)
