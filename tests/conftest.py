import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def fixtures_dir():
    return FIXTURES_DIR


@pytest.fixture
def run_ijump(tmp_path):
    """Run the ijump CLI as a subprocess against the tiny fixture set.

    Returns the completed process. Output/work directories are created fresh
    under tmp_path for each call, so runs don't interfere with each other or
    leave anything behind in the repo.
    """
    def _run(estimation_mode, extra_args=()):
        outdir = tmp_path / f"out_{estimation_mode}"
        workdir = tmp_path / f"wd_{estimation_mode}"
        args = [
            sys.executable,
            str(REPO_ROOT / "ijump.py"),
            "--aln", str(FIXTURES_DIR / "tiny.bam"),
            "--ref", str(FIXTURES_DIR / "tiny.fna"),
            "--gff", str(FIXTURES_DIR / "tiny.gff"),
            "--isel", str(FIXTURES_DIR / "is_coords.txt"),
            "--outdir", str(outdir),
            "--wd", str(workdir),
            "--estimation_mode", estimation_mode,
            *extra_args,
        ]
        result = subprocess.run(
            args,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        return result, outdir

    return _run
