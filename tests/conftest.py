import subprocess

import pytest
from golden_support import FIXTURES_DIR, REPO_ROOT


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

    def _run(estimation_mode, extra_args=(), isel="is_coords_clustered.txt"):
        """estimation_mode=None omits --estimation_mode entirely, exercising
        the CLI's default rather than an explicitly-passed value.

        ``isel`` names the IS table fixture to run against. The default carries
        a cluster column, which precise mode requires (isfinder-annotation 06);
        pass ``is_coords.txt`` for the legacy four-column table."""
        label = estimation_mode if estimation_mode is not None else "default"
        outdir = tmp_path / f"out_{label}"
        workdir = tmp_path / f"wd_{label}"
        args = [
            "ijump",
            "run",
            "--aln",
            str(FIXTURES_DIR / "tiny.bam"),
            "--ref",
            str(FIXTURES_DIR / "tiny.fna"),
            "--gff",
            str(FIXTURES_DIR / "tiny.gff"),
            "--isel",
            str(FIXTURES_DIR / isel),
            "--outdir",
            str(outdir),
            "--wd",
            str(workdir),
            *(["--estimation_mode", estimation_mode] if estimation_mode is not None else []),
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
