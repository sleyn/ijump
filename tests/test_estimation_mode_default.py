"""Regression test for ticket 14: omitting --estimation_mode must resolve to
a genuine EstimationMode.AVERAGE member, not the str(EstimationMode.AVERAGE)
literal "EstimationMode.AVERAGE" that argparse produced before the fix (see
the ticket for the full argparse-default-goes-through-type mechanics).

Before the fix, ISClipped.run()'s `mode == EstimationMode.AVERAGE` /
`mode == EstimationMode.PRECISE` dispatch matched neither branch when the
flag was omitted, so the pipeline silently no-op'd (exit 0, only reads.txt
written). test_empty_run_writes_full_file_set_with_no_rows below is the
end-to-end counterpart: it drives the CLI without --estimation_mode and
asserts the full average-mode file set is written.
"""
import pandas as pd
import pytest

from conftest import REPO_ROOT

from ijump.ijump import parse_args
from ijump.isclipped import EstimationMode


def test_omitted_estimation_mode_resolves_to_average_member():
    args = parse_args([])

    assert args.estimation_mode == EstimationMode.AVERAGE
    assert isinstance(args.estimation_mode, EstimationMode)


@pytest.mark.parametrize("estimation_mode", ["average", "precise"])
def test_explicit_estimation_mode_still_resolves_to_correct_member(estimation_mode):
    args = parse_args(["--estimation_mode", estimation_mode])

    assert args.estimation_mode == EstimationMode(estimation_mode)
    assert isinstance(args.estimation_mode, EstimationMode)


def test_omitted_estimation_mode_writes_full_average_mode_output_set(run_ijump):
    result, outdir = run_ijump(estimation_mode=None)

    assert result.returncode == 0, result.stderr

    for filename in ["reads.txt", "ijump_junctions.txt", "ijump_sum_by_reg.txt", "ijump_report_by_is_reg.txt"]:
        file_path = outdir / filename
        assert file_path.exists(), f"{filename} was not written"
        table = pd.read_csv(file_path, sep="\t")
        assert len(table) == 0, f"{filename} has data rows"
