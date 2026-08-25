# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Development uses `uv`, but `uv sync`/`uv run` currently cannot resolve this project's runtime deps
from PyPI alone (`pysamstats` hard-pins an unbuildable `pysam==0.15.4` — see README's "uv" section
for the full chain). Until ticket 04's conda packaging lands, use conda for anything that needs the
runtime deps installed, and `uv`/plain `pip` only for lint, tests, and build:

```
conda env create -f environment.yml && conda activate ijump
pip install -e . --no-deps        # editable install, resolves via conda's pysam/pysamstats
```

- Run all tests: `pytest` (from repo root; `pytest.ini` sets `testpaths = tests`)
- Run a single test file/case: `pytest tests/test_isclipped.py::test_name`
- Characterization goldens live in `tests/goldens/` (see its README). The end-to-end tier is
  marked `e2e`, needs the large inputs in `Test/` (or `$IJUMP_E2E_DATA`), and skips without
  them; deselect it with `pytest -m "not e2e"`. Re-pin an intended change with
  `python tests/regenerate_goldens.py` and review the resulting diff.
- Lint/format: `ruff check src/ijump tests` / `ruff format src/ijump tests` (mypy: `mypy src/ijump tests`)
- Install the pre-commit hook once per clone: `pre-commit install` (runs ruff --fix, ruff format, mypy on `git commit`); run over everything with `pre-commit run --all-files`
- `uv build` works standalone (doesn't need the runtime deps to resolve) for producing a wheel/sdist
- Structural search/codemods: `ast-grep scan .` runs this repo's own rules in `rules/`; `ast-grep test` checks them; see `docs/agents/ast-grep.md` for repo-specific recipes (state-coupling matrix, outline usage) before reaching for grep on a Python identifier

Lint/mypy in CI (`.github/workflows/lint.yml`) only cover `src/ijump/` and `tests/` — other
directories (`simulation/`, `rule-tests/`, root-level scripts) are out of scope and unclean.
CI does not run pytest (no working dependency-install recipe there yet); rely on local runs.

## Architecture

iJump detects Insertion Sequence (IS) element rearrangements in evolved bacterial populations from
short-read alignments, using **soft-clipped reads** as the signal: a read spanning an insertion
junction has one part mapped to the reference and the clipped remainder matching the IS element.
Full algorithm write-up: `docs/algorithm.md`. Domain vocabulary (avoid drifting to synonyms it
explicitly rejects): `CONTEXT.md`.

**Entry point**: `ijump` console script → `src/ijump/cli.py` dispatches to one of three subcommands,
each still parsing its own argv independently (this repo predates the unified CLI and the dispatch
layer is deliberately thin — it does not reinterpret any target's flags):
- `ijump run` → `ijump.ijump:main` — the main detection pipeline
- `ijump combine-results` → `ijump.combine_results:main` — merges per-sample report tables
- `ijump isfinder-db-parse` → `ijump.isfinder_db_parcer:main` — parses ISFinder BLAST outfmt-6
  output into the IS table (`is_table.py`), grouping copies of one element by sequence
  similarity (`is_clustering.py`)
- `ijump migrate-is-table` → `ijump.migrate_is_table:main` — annotates an existing
  (legacy four-column) IS table in place of regenerating it, preserving its coordinates
- `ijump isescan-convert` → `ijump.isescan_convert:main` — converts ISEScan's `.tsv`
  results into an IS table. iJump **reads** ISEScan output and never invokes ISEScan

The IS-table back-ends differ only in where the four locus columns come from; everything
after that — clustering and the origin-spanning flags — is `is_annotation.annotate_and_cluster`,
shared so they cannot drift on what a cluster is.

**Pipeline core**: `ISClipped` in `src/ijump/isclipped.py` drives both workflows and owns most
detection state across its methods (see `docs/agents/ast-grep.md`'s state-coupling matrix for which
attributes are single-write-then-read vs. mutated from multiple methods — relevant before extracting
anything from this class). It composes several single-responsibility modules rather than doing the
work inline:
- `clipped_read_search.py` — turns an alignment into a BLAST hit table of candidate junction
  positions (`direction=1`: IS→Ref forward search; `direction=0`: Ref→IS backward search, precise
  mode only)
- `junction_pairing.py` — pairs left/right junctions into IS element insertion positions (precise mode)
- `frequency_estimation.py` — estimates population frequency from paired junctions and depth
- `region_summary.py` — average mode's per-annotated-region report generation
- `circos.py` — renders finished detection results into Circos input files (`-c/--circos` flag);
  has no role in detection itself
- `gff.py` — a custom GFF reader/parser scoped to PATRIC/PROKKA-style GFFs specifically (see the
  GFF format requirements in README's Input section — other GFF flavors need reformatting)

**Two workflows**, selected by `--estimation_mode`:
- **Average** (`docs/Average.md`) — estimates IS frequency per gene/intergenic region from average
  region coverage vs. clipped-read count; lower accuracy, higher sensitivity; output centers on
  `ijump_report_by_is_reg.txt`
- **Precise** (`docs/Precise.md`) — localizes insertion coordinates via junction pairing before
  frequency estimation; more accurate, less tested; output centers on `ijump_junction_pairs.txt`

Both modes share the initial clipped-read collection and BLAST steps before diverging.

**Coverage accumulation** (`ISClipped.average_depth`): pure-`pysam` per-read CIGAR/span accumulator
(no `pysamstats`, no pileup). Excludes `UNMAP|SECONDARY|QCFAIL|DUP` reads (htslib's own default
pileup filter) but deliberately keeps `SUPPLEMENTARY` reads — iJump is a transposon-insertion caller
and supplementary alignments of clipped reads are exactly what it's looking for.

**Multi-sample comparison**: `combine_results.py` merges several samples' `ijump_report_by_is_reg.txt`
files (named `ijump_<sample>.txt`, all copied into one directory) into a single summary table, with
optional GFF-based functional annotation and clonal-sample comparison.

## Repo conventions

- Issues/specs are markdown files under `.scratch/<feature-slug>/`, not GitHub issues — see
  `docs/agents/issue-tracker.md`. Triage state uses the label vocabulary in `docs/agents/triage-labels.md`.
- Domain terms belong in `CONTEXT.md`; check it before introducing new vocabulary for an existing concept.
