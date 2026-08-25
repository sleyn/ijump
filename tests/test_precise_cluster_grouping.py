"""Precise mode groups junctions by the IS table's cluster column.

Which junctions can be paired into one insertion is a question about which loci
a clipped read cannot tell apart, and that is what the cluster column records
(isfinder-annotation 06). It used to be inferred by stripping a trailing ``_1``
off the IS name -- a suffix the parser invented in the first place, which gets
the reference genome wrong in both directions: the two ``IS17`` fragments do not
align to each other at all, while their actual parent ``ISAba12_1`` was held
apart from both.
"""

import pandas as pd
import pytest
from fake_alignment import FakeAlignment

from ijump import is_table
from ijump.isclipped import EstimationMode, ISClipped

# One insertion of an ISAba12-like element: a left junction whose clipped part
# matched the ``IS17_1`` fragment and a right junction 5 bp away that matched
# ``ISAba12_1``. The two are one element, so the two junctions are one pair.
JUNCTIONS = pd.DataFrame(
    {
        "IS name": ["IS17_1", "ISAba12_1", "ISAba1_1"],
        "Chrom": ["contig_1"] * 3,
        "Position": [4000, 4005, 7000],
        "Orientation": ["left", "right", "left"],
    }
)

IS_TABLE = pd.DataFrame(
    [
        ["IS17_1", "contig_1", "100", "244", "IS5", "IS903", "ISAba12", "98.6", "no", ""],
        ["ISAba12_1", "contig_1", "1000", "2039", "IS5", "IS903", "ISAba12", "98.5", "no", ""],
        ["ISAba1_1", "contig_1", "3000", "4179", "IS4", "IS10", "ISAba1", "100", "no", ""],
    ],
    columns=list(is_table.COLUMNS),
)


def _pipeline(is_table_frame, tmp_path):
    isc = ISClipped(FakeAlignment(), "ref", "unused.gff", str(tmp_path))
    isc.is_table = is_table_frame
    isc.junctions = JUNCTIONS.copy()
    return isc


def test_junctions_of_one_cluster_pair_across_is_names(tmp_path):
    isc = _pipeline(IS_TABLE, tmp_path)

    isc.search_insert_pos()

    aba12 = isc.pairs_df.query('IS_name == "ISAba12"')
    assert len(aba12) == 1
    assert aba12.iloc[0]["Position_l"] == 4000
    assert aba12.iloc[0]["Position_r"] == 4005


def test_pairs_are_labelled_with_the_cluster_not_the_stripped_name(tmp_path):
    isc = _pipeline(IS_TABLE, tmp_path)

    isc.search_insert_pos()

    assert set(isc.pairs_df["IS_name"]) == {"ISAba12", "ISAba1"}


def test_precise_run_on_a_table_without_clusters_fails_before_any_search(tmp_path):
    """The error is raised up front, not after a couple of minutes of BLAST."""
    legacy = IS_TABLE.copy()
    legacy["cluster"] = ""
    isc = _pipeline(legacy, tmp_path)

    with pytest.raises(is_table.MissingClusterColumn) as excinfo:
        isc.run(EstimationMode.PRECISE)

    assert is_table.MIGRATE_SUBCOMMAND in str(excinfo.value)


def test_cli_precise_run_on_a_legacy_table_reports_the_remedy(run_ijump):
    """Through the CLI the failure is a message, not a traceback: a table with
    no cluster column is a setup problem with a documented remedy."""
    result, _ = run_ijump("precise", isel="is_coords.txt")

    assert result.returncode == 1
    assert is_table.MIGRATE_SUBCOMMAND in result.stdout


def test_cli_average_run_still_accepts_a_legacy_table(run_ijump):
    """Only precise mode groups on the cluster, so average mode is not
    retired here -- that is isfinder-annotation 07's call to make."""
    result, _ = run_ijump("average", isel="is_coords.txt")

    assert result.returncode == 0, result.stderr
