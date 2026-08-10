"""CLI-level test for ticket 05: an invalid --estimation_mode value must be
rejected by argparse at parse time, not silently accepted and left to fall
through every EstimationMode comparison in the workflow.

Unlike run_ijump's fixture runs, no other flags are needed here: argparse
validates choices for an option as soon as it's parsed, before checking
that other (unrelated) arguments were supplied at all.
"""
import subprocess
import sys

from conftest import REPO_ROOT


def test_invalid_estimation_mode_rejected_at_parse_time():
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "ijump.py"), "--estimation_mode", "bogus"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "invalid choice" in result.stderr
