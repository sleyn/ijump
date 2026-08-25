"""``isfinder-db-parse`` flags the copy the assembler broke (isfinder-annotation 05).

The detection rule is pinned in ``test_origin_spanning.py``; this covers the
wiring -- that the parser reads the contig lengths out of the reference it was
given and puts the answer in the table. ``NODE_2`` of
``Test/A_baumannii_assembly.fna`` carries the case: 148137 bases with an
ISAba12-like copy across the seam.
"""

import golden_support
import pandas as pd
import pytest


@pytest.fixture(autouse=True)
def requires_blast():
    reasons = golden_support.missing_parser_requirements()
    if reasons:
        pytest.skip("; ".join(reasons))


@pytest.fixture
def table(tmp_path):
    result = golden_support.run_isfinder_db_parse(tmp_path)
    assert result.returncode == 0, result.stderr
    return pd.read_csv(
        tmp_path / golden_support.ISFINDER_TABLE_NAME, sep="\t", dtype=str, keep_default_na=False
    )


def row(table, is_name):
    (found,) = table[table["is_name"] == is_name].to_dict("records")
    return found


def test_the_two_is17_fragments_are_flagged_and_share_an_element_id(table):
    first, second = row(table, "IS17_1"), row(table, "IS17_2")

    assert first["wraps_origin"] == "yes"
    assert second["wraps_origin"] == "yes"
    assert first["element_id"] == second["element_id"] != ""


def test_nothing_else_on_the_genome_is_flagged(table):
    """Every other locus sits inside its contig; ``ISAlw13_1`` is 553 bases from
    the start of NODE_3 with nothing at the far end, which is not a break."""
    flagged = set(table.loc[table["wraps_origin"] == "yes", "is_name"])

    assert flagged == {"IS17_1", "IS17_2"}
    assert set(table.loc[table["wraps_origin"] == "no", "element_id"]) == {""}


def test_both_coordinate_rows_survive_with_start_before_stop(table):
    """The halves are never joined into one ``start > stop`` row: too much downstream
    assumes a row's start precedes its stop."""
    assert (table["start"].astype(int) <= table["stop"].astype(int)).all()

    first, second = row(table, "IS17_1"), row(table, "IS17_2")
    assert (first["start"], first["stop"]) == ("147994", "148137")
    assert (second["start"], second["stop"]) == ("2", "77")
