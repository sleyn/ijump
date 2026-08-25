"""End-to-end test for the empty-run file-set invariant (ticket 02): a run
that finds nothing writes the same files as a run that finds something,
with headers and zero data rows, in both estimation modes.

With the tiny.bam fixture (reads that align but are never soft-clipped),
execution raises NoInsertionsFound at check_data_presence_in_df
(ijump.py:35, called at ijump.py:193), before runblast (ijump.py:198) ever
runs. Because tiny.nsq is committed, check_blast_ref (ijump.py:53)
short-circuits and never shells out to makeblastdb, so this test needs no
BLAST+ install.

Before ticket 02, this exit path wrote only ijump_junction_pairs.txt and
skipped the downstream report files, silently dropping the sample out of
combine_results.py (see the ticket's "Why"). See git history for the prior
version of this test, which pinned that behaviour.
"""

import pandas as pd
import pytest

from ijump import report_provenance

FILES_BY_MODE = {
    "average": [
        "reads.txt",
        "ijump_junctions.txt",
        "ijump_sum_by_reg.txt",
        "ijump_report_by_is_reg.txt",
    ],
    "precise": [
        "reads.txt",
        "ijump_junctions.txt",
        "ijump_sum_by_reg.txt",
        "ijump_report_by_is_reg.txt",
        "ijump_junction_pairs.txt",
    ],
}


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_empty_run_exits_cleanly(run_ijump, estimation_mode):
    result, outdir = run_ijump(estimation_mode)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_empty_run_writes_full_file_set_with_no_rows(run_ijump, estimation_mode):
    result, outdir = run_ijump(estimation_mode)
    assert result.returncode == 0, result.stderr

    for filename in FILES_BY_MODE[estimation_mode]:
        file_path = outdir / filename
        assert file_path.exists(), f"{filename} was not written"
        table = _read(file_path)
        assert len(table) == 0, f"{filename} has data rows"


# The two files combine_results merges carry a leading line naming the IS table
# the run was annotated against (isfinder-annotation 07); the rest are plain TSVs.
STAMPED = {"ijump_report_by_is_reg.txt", "ijump_junction_pairs.txt"}


def _read(file_path):
    if file_path.name in STAMPED:
        table, _ = report_provenance.read_report(file_path)
        return table
    return pd.read_csv(file_path, sep="\t")
