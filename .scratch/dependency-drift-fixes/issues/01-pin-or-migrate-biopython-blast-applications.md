# 01 — Stop clipped_read_search.py from breaking on current Biopython

Status: done

**What to build:** `src/ijump/clipped_read_search.py:21` does
`from Bio.Blast.Applications import NcbiblastnCommandline`. Biopython
deprecated `Bio.Blast.Applications` around 1.79/1.80 and removed it entirely
in later releases; `pyproject.toml` pins plain `"biopython"` with no version
ceiling, so a fresh `uv`/`pip` resolve grabs the latest release and this
import fails at collection time — taking down every test module that
transitively imports `isclipped.py` (i.e. most of the suite), not just BLAST
tests. Fix so the BLAST call works against whatever Biopython a fresh
install resolves to, without needing a hand-picked `biopython<1.80` in a
throwaway venv to get `pytest` running at all.

Two acceptable approaches — pick one:

1. Pin `biopython<1.80` (or another ceiling that still has
   `Bio.Blast.Applications`) in `pyproject.toml`'s `dependencies`, matching
   whatever `environment.yml`/conda already resolves to for this package if
   it's already pinned there.
2. Migrate the `NcbiblastnCommandline` call to a direct
   `subprocess.run(["blastn", ...])` invocation — Biopython's own documented
   replacement for the removed `Applications` wrappers — so newer Biopython
   (and the eventual removal of the deprecated module from PyPI's latest)
   doesn't matter.

**Blocked by:** None — can start immediately

- [x] `uv venv` + `uv pip install -e . pytest` (no biopython version pin of
      your own) resolves and installs cleanly on a supported Python version
- [x] `pytest` collects and runs without `ModuleNotFoundError: No module
      named 'Bio.Blast.Applications'` in that fresh environment
- [x] The BLAST invocation itself (`clipped_read_search.search`) still
      produces the same output as before — existing BLAST-dependent tests
      (`tests/test_clipped_read_search.py`, `tests/test_check_blast_ref.py`)
      pass unchanged
- [x] `ruff check`/`mypy` on `src/ijump` and `tests` still pass

## Comments

Took approach 2: migrated `run_blast_subprocess` off `NcbiblastnCommandline`
onto a direct `subprocess.run(["blastn", ...], check=True)` call — Biopython's
own documented replacement for the removed `Applications` wrappers. Biopython
is no longer imported anywhere in `src/ijump`, so it was dropped entirely from
`pyproject.toml`'s `dependencies`, `environment.yml`, and README's dependency
list rather than left as dead weight.

Verified against `/private/tmp/ijump-venv2` (pysam 0.24.0, pytest 8.4.2, blastn
on PATH): `tests/test_clipped_read_search.py` and `tests/test_check_blast_ref.py`
pass both with and without biopython installed at all. Full suite passes except
4 pre-existing `test_average_depth_pysamstats_parity.py` failures (that venv has
no `pysamstats` installed — unrelated to this ticket; `pysamstats` was already
dropped from the runtime deps per ticket 16).

Post-review touch-up: added the same `FileNotFoundError`/`CalledProcessError`
handling `ijump.py`'s `check_blast_ref`/`makeblastdb_command` already uses for
its own BLAST+ subprocess call, so a missing `blastn` binary or a failing run
now logs the same actionable message instead of a bare traceback.
