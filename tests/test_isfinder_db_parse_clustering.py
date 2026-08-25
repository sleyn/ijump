"""``isfinder-db-parse`` end to end over the committed fixtures (isfinder-annotation 04).

The clustering rule itself is pinned in ``test_is_clustering.py`` against the
measured alignments; these tests cover the wiring around it -- the reference
input, the threshold flags, and the invariant the column has to hold.

Inputs are the ISFinder BLAST output of ``Test/A_baumannii_assembly.fna`` and a
masked stand-in for that assembly (see
``tests/fixtures/isfinder/make_reference_fixture.py``), both a few kilobytes.
BLAST+ has to be on PATH.
"""

import golden_support
import pandas as pd
import pytest


@pytest.fixture(autouse=True)
def requires_blast():
    reasons = golden_support.missing_parser_requirements()
    if reasons:
        pytest.skip("; ".join(reasons))


def parse(tmp_path, extra_args=()):
    result = golden_support.run_isfinder_db_parse(tmp_path, extra_args=extra_args)
    assert result.returncode == 0, result.stderr
    return pd.read_csv(tmp_path / golden_support.ISFINDER_TABLE_NAME, sep="\t", dtype=str)


def clusters(table):
    return dict(zip(table["is_name"], table["cluster"]))


def test_cluster_names_are_unique_across_the_table(tmp_path):
    """``ISClipped.is_coords`` is a dict keyed on the grouping name and precise
    mode groups on it, so a duplicate would silently swallow a locus rather than
    just look untidy."""
    table = parse(tmp_path)

    names = set(table["cluster"])
    assert len(names) == table["cluster"].nunique(dropna=False)
    assert "" not in names and not table["cluster"].isna().any()


def test_every_row_gets_a_cluster(tmp_path):
    table = parse(tmp_path)

    assert len(table) == len(table.dropna(subset=["cluster"]))


def test_the_is17_fragments_and_isaba12_share_a_cluster(tmp_path):
    """The three rows are one ISAba12-like copy plus its own two contig-edge
    fragments; ``ISAba53_1`` is 83% away and a different element."""
    assigned = clusters(parse(tmp_path))

    assert assigned["IS17_1"] == assigned["ISAba12_1"]
    assert assigned["IS17_2"] == assigned["ISAba12_1"]
    assert assigned["ISAba12_1"] == "ISAba12"
    assert assigned["ISAba53_1"] != assigned["ISAba12_1"]


def test_the_three_isaba11_copies_share_a_cluster(tmp_path):
    assigned = clusters(parse(tmp_path))

    assert assigned["ISAba11_1"] == "ISAba11"
    assert assigned["ISAba11_2"] == "ISAba11"
    assert assigned["ISAba11_3"] == "ISAba11"


def test_lowering_the_identity_flag_merges_isaba53_into_isaba12(tmp_path):
    """The 83% pair is the flag's own test case: it is deliberately out at the
    default and in once the floor drops below it."""
    assigned = clusters(parse(tmp_path, extra_args=["--cluster-identity", "80"]))

    assert assigned["ISAba53_1"] == assigned["ISAba12_1"]


def test_the_coverage_flag_reaches_the_linkage_rule(tmp_path):
    """No pair on this genome aligns over only *part* of the shorter locus --
    every one of them covers it whole -- so nothing short of an unreachable
    coverage moves the answer here. That still pins the flag as wired: raise it
    past 1 and every locus stands alone.

    The threshold's actual behaviour is pinned on measured numbers in
    ``test_is_clustering.py``.
    """
    assigned = clusters(parse(tmp_path, extra_args=["--cluster-coverage", "1.01"]))

    assert len(set(assigned.values())) == len(assigned)


def test_a_reference_that_does_not_carry_the_called_contigs_is_an_error(tmp_path):
    """Pointing at the wrong genome would otherwise cluster whatever sequence
    happened to sit at those coordinates."""
    wrong = tmp_path / "wrong.fna"
    wrong.write_text(">some_other_contig\n" + "ACGT" * 100 + "\n")

    result = golden_support.run_isfinder_db_parse(tmp_path, reference=wrong)

    assert result.returncode != 0
    assert "is not in" in result.stderr
