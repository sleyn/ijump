# 01 - Move to a `src/ijump` package layout

## Context

The repo has 11 flat top-level `.py` modules (`circos.py`, `clipped_read_search.py`,
`combine_results.py`, `frequency_estimation.py`, `gff.py`, `ijump.py`, `isclipped.py`,
`isfinder_db_parcer.py`, `isfinder_parser.py`, `junction_pairing.py`, `region_summary.py`)
that import each other by bare name, relying on pytest's implicit rootdir-on-sys.path.
There is no `pyproject.toml`.

## Checklist

- [x] Create `src/ijump/`, move all 11 modules into it, add `src/ijump/__init__.py`.
- [x] Rewrite every intra-package import to be package-relative (`from ijump import x` /
      `from ijump.x import y`). Known edges: `combine_results.py` (`gff`), `ijump.py`
      (`isclipped`, `circos`), `isclipped.py` (`gff`, `frequency_estimation`,
      `clipped_read_search`, `region_summary`, `junction_pairing`). Re-grepped at
      implementation time; no other intra-package edges were found.
      `tests/test_average_depth_pysamstats_vs_count_coverage.py` (added by ticket 16)
      imported `from isclipped import ISClipped` — updated to `from ijump.isclipped
      import ISClipped`.
- [x] Wrap `isfinder_parser.py` and `isfinder_db_parcer.py`'s top-level
      argument-parsing code in a `main()` function behind
      `if __name__ == "__main__":`, so importing them as package submodules
      doesn't try to parse `sys.argv`.
- [x] Add `pyproject.toml` (PEP 621 `[project]` table): package name `ijump`,
      src layout, a build backend, and Python dependencies implied by imports.
      Did NOT declare BLAST+ (non-Python system binary, handled gracefully at
      runtime by `ijump.py`'s `check_blast_ref`).
- [x] No console-script entry points added (out of scope for this ticket).
      Package is invocable as `python -m ijump.ijump ...` (and
      `python -m ijump.<module>` for the other CLI scripts).
- [x] Updated all 14 test files under `tests/` to import via
      `from ijump.<module> import ...` / `import ijump.<module> as <module>`
      instead of bare names. Updated `tests/conftest.py`'s `run_ijump` fixture
      and `tests/test_estimation_mode_validation.py`'s direct subprocess call
      to invoke `python -m ijump.ijump` instead of a hardcoded script path.
- [x] `pytest` resolves `ijump.*` imports without manual `PYTHONPATH` fiddling,
      via `pip install -e .`. Documented in the README's new "Development setup"
      section under Installation.
- [x] Left non-Python directories/files at the repo root alone: `Test/`,
      `Example files/`, `ijump_data/`, `img/`, `simulation/`, `circos.conf`.

## Verification

- `pip install -e .` succeeds from a clean state (uninstall + reinstall
  exercised, in the conda env at `/Users/sleyn/miniconda3/envs/bioinfo`).
- `pytest` (same conda env's Python): 45 passed, 1 failed
  (`test_read_count_mtx_rejects_invalid_orientation` — pre-existing,
  references a since-removed method, not in scope for this ticket).
- No bare top-level `.py` modules remain at the repo root.

## Comments

- **Build backend:** chose `setuptools` (`setuptools>=61` for PEP 621
  `[project]` table support) over `hatchling` — it's the more conservative,
  widely-available choice and the target conda env's pip could install a
  fresh isolated-build copy of it without any extra configuration.
- **`requires-python`:** set to `>=3.7` to match the conda env
  (`bioinfo`, Python 3.7.6) actually used to run the test suite; `pysam`/
  `pysamstats` are conda-only on most platforms, so that env's Python version
  effectively pins this.
- **Dependencies:** the ticket text named only `pandas`, `pysam`,
  `pysamstats`, `biopython` as implied by imports, but a re-grep at
  implementation time found `isclipped.py` also directly imports `numpy`
  and `from sklearn.cluster import AgglomerativeClustering`. Both were
  added to `[project.dependencies]` (as `numpy` and `scikit-learn`) since
  omitting them would make `pip install -e .` succeed but leave the package
  unimportable in a genuinely clean environment.
- **pythonpath mechanism:** picked `pip install -e .` (editable install)
  over a `pytest.ini` `pythonpath` setting, and documented it in a new
  README "Development setup" subsection. Note: the conda env's pip (20.2.4)
  was too old to support PEP 660 editable installs from a `pyproject.toml`
  with no `setup.py`; it was upgraded in-place to pip 24.0 (last version
  supporting Python 3.7) to make `pip install -e .` work.
- **Script invocation gotcha:** running a submodule directly by path (e.g.
  `python src/ijump/ijump.py` or `python src/ijump/combine_results.py`) puts
  `src/ijump/` first on `sys.path`. Since `ijump.py` lives in that same
  directory, `from ijump import ...` / `from ijump.x import y` then resolves
  `ijump` to that single file instead of the installed package, raising
  `ModuleNotFoundError: No module named 'ijump.isclipped'; 'ijump' is not a
  package`. Fixed by invoking everything as `python -m ijump.<module>`
  instead (per ticket item 5's suggested invocation form) — updated
  `tests/conftest.py`, `tests/test_estimation_mode_validation.py`, and all
  README usage examples accordingly.
- **`circos.conf` lookup:** `circos.py`'s `write_files` locates `circos.conf`
  via `os.path.dirname(os.path.realpath(__file__))`, which used to be the
  repo root when `circos.py` lived there. After the move, `circos.py` lives
  two directories below the repo root (`src/ijump/circos.py`), so this was
  updated to walk up two more directories to keep finding the
  repo-root `circos.conf` (left in place per the ticket's out-of-scope
  list). This is a pure path fix required by the move, not a behavior
  change — `tests/test_circos.py` now passes again with it.
- Also added `*.egg-info/`, `build/`, `dist/` to `.gitignore` for the new
  packaging build artifacts.
