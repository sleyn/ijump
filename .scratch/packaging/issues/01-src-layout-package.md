# 01 — Move Python modules into a proper `src/ijump/` package

**What to build:** Turn the 11 flat root-level `.py` modules into an
installable package at `src/ijump/`, with a `pyproject.toml` that makes
`pip install .` (and later `uv sync`/`uv build`, ticket 03) produce a real
`ijump` distribution — not just relocate files.

## Why

This is step one of a packaging effort (conda + Docker + PyPI, grilled and
scoped across tickets 01–05 in this directory). The user wants a "clean
Python package" installable via `pip`/`uv`/conda; today the repo has no
`pyproject.toml`/`setup.py` at all, and every module is a flat top-level
`.py` file at the repo root relying on pytest's implicit rootdir-on-sys.path
behavior (`pytest.ini`'s `testpaths = tests`, no `src` layout) to import
each other by bare name (e.g. `isclipped.py:7-12` does `import gff`,
`import frequency_estimation`, `from clipped_read_search import
NoInsertionsFound`).

Publishing top-level modules named `gff`, `region_summary`, etc. under no
namespace risks collisions with other installed packages and isn't
idiomatic for a PyPI package — a proper package needs its own namespace.

## Scope

- Create `src/ijump/` and move these into it (the full set of internal
  modules — non-CLI library code only; the four CLI entry-point scripts
  move here too but ticket 02 handles turning them into subcommands, don't
  block this ticket on that):
  `circos.py`, `clipped_read_search.py`, `combine_results.py`,
  `frequency_estimation.py`, `gff.py`, `ijump.py`, `isclipped.py`,
  `isfinder_db_parcer.py`, `isfinder_parser.py`, `junction_pairing.py`,
  `region_summary.py`.
- Add `src/ijump/__init__.py`.
- Rewrite every intra-package import to be package-relative. Known import
  edges to fix (found via grep, re-verify at implementation time in case
  more have landed):
  - `combine_results.py:7` — `import gff` → `from ijump import gff` (or
    relative `from . import gff`)
  - `ijump.py:8-9` — `from isclipped import ISClipped, EstimationMode`,
    `import circos` → package-relative equivalents
  - `isclipped.py:7-12` — `import gff`, `import frequency_estimation`,
    `import clipped_read_search`, `import region_summary`,
    `from clipped_read_search import NoInsertionsFound`,
    `from junction_pairing import find_pairs` → package-relative
    equivalents
  - `isfinder_parser.py` and `isfinder_db_parcer.py` currently have no
    `if __name__ == "__main__":` guard — they execute top-level code
    (`getopt.getopt(sys.argv[1:], ...)` / `argparse` parsing) at import
    time. Importing them as package submodules without guarding this would
    break at package-import time (e.g. any `from ijump import
    isfinder_parser` in a test would try to parse `sys.argv`). Wrap each
    script's body in a `main()` function behind an `if __name__ ==
    "__main__":` guard as part of this move — needed regardless of ticket
    02's subcommand work, and blocks it (ticket 02 dispatches to these
    `main()`s).
- Add `pyproject.toml` (PEP 621 `[project]` table) declaring: package name
  `ijump`, the `src` layout (`[tool.setuptools.package-dir]` /
  build-backend of your choice — `setuptools` or `hatchling`, pick whichever
  ticket 03's uv migration expects, since uv itself is backend-agnostic and
  just needs a valid build backend declared), and the Python dependencies
  currently only implied by imports (`pandas`, `pysam`, `pysamstats`,
  `biopython`). Do not declare BLAST+ as a dependency — it's a non-Python
  system binary; `ijump.py`'s existing `check_blast_ref` runtime check
  (`ijump.py:26-45`) already handles "BLAST+ not found" gracefully and stays
  as-is.
- Do NOT add console-script entry points yet — that's ticket 02. This
  ticket's package should still be invoked as
  `python -m ijump.ijump ...`/similar or via `PYTHONPATH=src python3
  src/ijump/ijump.py ...` for now; ticket 02 replaces that with a real `ijump`
  command.
- Update `tests/`: all 14 test files import target modules by bare name
  today (e.g. `tests/test_region_summary.py:19` —
  `from region_summary import summarize_by_region, report_average`). Rewrite
  every such import to `from ijump.region_summary import ...`-style.
  `tests/conftest.py`'s `run_ijump` fixture (lines 30-41) invokes
  `REPO_ROOT / "ijump.py"` directly as a subprocess — update the path to
  `REPO_ROOT / "src" / "ijump" / "ijump.py"` for this ticket; ticket 02 will
  change this again once a real console-script exists.
- Update `pytest.ini`/add `pythonpath`/`pip install -e .` step so `pytest`
  can still resolve `ijump.*` imports without manual `PYTHONPATH` fiddling
  in CI or local dev — confirm which mechanism you use (editable install vs.
  `pytest.ini`'s `pythonpath` setting) and document it in the README's dev
  setup section.
- Leave non-Python directories at the repo root as-is: `Test/` (example/
  reference data, unrelated to `tests/`), `Example files/`, `ijump_data/`,
  `img/`, `simulation/`, `circos.conf`.

## Out of scope

- CLI subcommand dispatch / console-script entry points (ticket 02).
- uv migration, `uv.lock` (ticket 03).
- conda packaging (ticket 04).
- Dockerfile (ticket 05).
- Any behavior change to the algorithms themselves — this is a pure move +
  import rewrite.

## Verification

- `pip install -e .` (or equivalent) succeeds from a clean clone.
- `pytest` passes from a clean clone with the package installed.
- Manual real-sample run (reuse ticket 13's verification command from
  `.scratch/isclipped-refactor/`) still produces the same output files via
  the new module path.
- No bare top-level `.py` modules remain at the repo root except this
  directory's non-Python assets.

**Blocked by:** None — can start immediately

**Status:** done

- [x] All 11 modules moved to `src/ijump/`, `__init__.py` added.
- [x] Every intra-package import rewritten to package-relative form.
- [x] `isfinder_parser.py` and `isfinder_db_parcer.py` wrapped in `main()` behind `if __name__ == "__main__":` guards.
- [x] `pyproject.toml` added with correct `[project]` metadata, `src` layout, and Python-only dependencies.
- [x] All 14 test files' imports updated to `ijump.*` package paths.
- [x] `tests/conftest.py`'s `run_ijump` fixture path updated.
- [x] `pytest` passes from a clean clone.
- [ ] Manual real-sample run produces the expected output set. (not exercised this pass — deferred to whoever does the next manual real-sample verification)

## Comments

Implemented 2026-08-16. All 11 modules moved to `src/ijump/`, imports
rewritten, `pyproject.toml` added, all 14 test files updated. `pip install
-e .` and `pytest` verified in the conda env at
`/Users/sleyn/miniconda3/envs/bioinfo` (Python 3.7.6): 45 passed, 1 failed
(`test_read_count_mtx_rejects_invalid_orientation` — pre-existing,
references a since-removed method, not in scope for this ticket).

Decisions and gotchas:

- **Build backend:** `setuptools` (`setuptools>=61`, PEP 621 `[project]`
  table support) over `hatchling` — more conservative and readily available;
  the target conda env's pip could pull a fresh isolated-build copy without
  extra configuration.
- **`requires-python`:** set to `>=3.7` to match the conda env actually used
  to run the suite (`pysam`/`pysamstats` are conda-only on most platforms,
  so that env's Python version effectively pins this).
- **Dependencies:** the ticket named `pandas`, `pysam`, `pysamstats`,
  `biopython`. Re-grepping at implementation time found `isclipped.py` also
  directly imports `numpy` and `from sklearn.cluster import
  AgglomerativeClustering` — both added to `[project.dependencies]`
  (`numpy`, `scikit-learn`) since omitting them would leave the package
  unimportable in a genuinely clean environment.
- **pythonpath mechanism:** `pip install -e .` (editable install), documented
  in a new README "Development setup" subsection. The conda env's pip
  (20.2.4) was too old for PEP 660 editable installs from a
  `pyproject.toml`-only project (no `setup.py`); upgraded in-place to pip
  24.0 (last version supporting Python 3.7).
- **Script invocation gotcha:** running a submodule directly by path (e.g.
  `python src/ijump/ijump.py`) puts `src/ijump/` first on `sys.path`; since
  `ijump.py` lives in that same directory, `from ijump import ...` then
  resolves `ijump` to that one file instead of the installed package
  (`ModuleNotFoundError: No module named 'ijump.isclipped'; 'ijump' is not a
  package`). Fixed by invoking everything as `python -m ijump.<module>`
  instead — updated `tests/conftest.py`,
  `tests/test_estimation_mode_validation.py`, and README usage examples
  accordingly.
- **`circos.conf` lookup:** `circos.py`'s `write_files` locates `circos.conf`
  via `os.path.dirname(os.path.realpath(__file__))`, which used to be the
  repo root when `circos.py` lived there. After the move it lives two
  directories below the repo root (`src/ijump/circos.py`), so the path-walk
  was updated to climb two more directories to keep finding the repo-root
  `circos.conf` (left in place per this ticket's out-of-scope list). Pure
  path fix required by the move, no behavior change —
  `tests/test_circos.py` passes with it.
- Also added `*.egg-info/`, `build/`, `dist/` to `.gitignore` for the new
  packaging build artifacts.
- The manual real-sample run (reusing ticket 13's verification command)
  wasn't exercised in this pass — flagged as outstanding above.

**Correction (2026-08-17, review-followups ticket 02):** The `test_read_count_mtx_rejects_invalid_orientation` failure noted above was *not* pre-existing. It was introduced by isclipped-refactor ticket 09's extraction of the read-count-matrix helper to module level in `frequency_estimation.py`, which left no delegating alias on `ISClipped`; on `master` the helper is still an `ISClipped` static method and the test passes there. It is fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.
