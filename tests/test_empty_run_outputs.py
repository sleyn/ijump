"""Characterization test: pins what iJump does today when no clipped reads
are found, in both estimation modes. This is not a spec of correct
behaviour — see .scratch/isclipped-refactor/issues/02-*.md, which inverts the
three "does not exist" assertions below once the empty-clipped-reads exit
path is fixed to still produce the downstream report files.

With the tiny.bam fixture (reads that align but are never soft-clipped),
execution reaches check_data_presence_in_df (ijump.py:28, called at
ijump.py:189) and exits 0 after writing only ijump_junction_pairs.txt via
empty_pairs_out() (ijump.py:22), before runblast (ijump.py:194) ever runs.
Because tiny.nsq is committed, check_blast_ref (ijump.py:50) short-circuits
and never shells out to makeblastdb, so this test needs no BLAST+ install.
"""
import pandas as pd
import pytest


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_empty_run_exits_cleanly(run_ijump, estimation_mode):
    result, outdir = run_ijump(estimation_mode)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_empty_run_writes_only_pairs_file_with_no_rows(run_ijump, estimation_mode):
    result, outdir = run_ijump(estimation_mode)
    assert result.returncode == 0, result.stderr

    pairs_file = outdir / "ijump_junction_pairs.txt"
    assert pairs_file.exists()
    pairs_df = pd.read_csv(pairs_file, sep="\t")
    assert len(pairs_df) == 0


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_empty_run_does_not_write_downstream_reports(run_ijump, estimation_mode):
    """Pins current, INCORRECT behaviour: these files are silently skipped
    instead of being written empty. Ticket 02 inverts these assertions."""
    result, outdir = run_ijump(estimation_mode)
    assert result.returncode == 0, result.stderr

    assert not (outdir / "ijump_report_by_is_reg.txt").exists()
    assert not (outdir / "ijump_sum_by_reg.txt").exists()
    assert not (outdir / "ijump_junctions.txt").exists()
