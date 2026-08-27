# 01 — Remove the Circos output module entirely

**What to build:** iJump no longer offers Circos output. The `-c/--circos` flag,
`circos.py`, and every reference to Circos in docs and CLAUDE.md are gone. A user
running `ijump run --help` sees no trace of the feature ever having existed.

**Why:** Circos rendering has no role in detection itself (it consumes finished
detection results and writes Circos input files) and is effectively unused in
practice. Dropping it removes one of the "clump" call sites review-followups/09
was weighing, and is a clean breaking change to land in the 2.0.0 major bump
alongside the other breaking changes already in this release
(`.scratch/release-2.0.0/issues/01-cut-2-0-0-release.md`).

**Blocked by:** None — can start immediately.

- [x] `src/ijump/circos.py` deleted
- [x] `-c/--circos` CLI flag removed from `ijump run`'s argument parser
- [x] Any call sites invoking the Circos writer (from `ISClipped`/the pipeline
      orchestration) removed
- [x] `docs/algorithm.md`, README, and `CLAUDE.md`'s architecture section no
      longer mention Circos or `circos.py` (`docs/algorithm.md` never did;
      `README.md`, `CLAUDE.md`, `docs/Average.md` updated)
- [x] No test references `circos.py` or the `-c/--circos` flag; any Circos-only
      tests are deleted rather than left failing
- [x] `ruff check`/`ruff format --check`/`mypy` clean on `src/ijump/`
- [x] Full `pytest` run has no new failures relative to the pre-removal baseline
- [x] A CHANGELOG entry records the removal as a breaking change (2.0.0)

## Comments
