# 12 — Close the private-method leak between `ISClipped` and `circos.write_files`

Status: ready-for-agent
Blocked by: —

## Why

`ijump.py:150` calls `circos.write_files(..., is_processing._av_depth, ...)` — reaching into
`ISClipped`'s leading-underscore method (`isclipped.py:667-673`) and handing it across the
file seam as a callback. Every other parameter `circos.write_files` takes is already plain
data (`report_table`, `sum_by_region`, `is_coords`, `ref_len`, `gff_ann_pos`); `_av_depth`
alone is private implementation reaching across a module boundary — the interface `circos.py`
actually depends on isn't what its signature says, it's whatever `ISClipped` happens to expose.

Architecture review candidate 02 (`/improve-codebase-architecture`, 2026-08-13), grilled to a
shared understanding the same session. Vocabulary: `/codebase-design`.

Decisions from grilling, all load-bearing for the scope below:

1. **Interface shape stays a lazy callback**, not a precomputed dict. `circos.write_files`
   calls `av_depth` once per `sum_by_region` row and once per `gff_ann_pos` region — many
   distinct `(chrom, start, stop)` triples computed via `pysamstats.load_coverage` against the
   live BAM. `report_average` (`isclipped.py:696`) already calls the same method. Precomputing
   would require `ISClipped` to know circos's region set upfront and re-couple the two modules;
   the lazy callback keeps `@lru_cache` doing the work of avoiding duplicate `pysamstats` calls
   across both callers.
2. **Name**: `average_depth` — drops the leading underscore, keeps the existing name
   otherwise. Matches the domain vocabulary already in use (`report_table`'s `Depth` column).
3. **Scope is narrow**: no `DetectionResult` value object. That's architecture review
   candidate 01, a separate, larger piece of work — not a blocker for this ticket.
4. **`@lru_cache` is untouched.** `ISClipped` is a single per-run object discarded after
   `main()` returns; the cache is exactly what avoids recomputing `pysamstats` for regions
   `report_average` and `circos.write_files` both query. Not part of this ticket's problem.
5. **`circos.py` and `tests/test_circos.py` need no change.** `circos.write_files`'s
   parameter is already named `av_depth` and treated as an opaque callable; `test_circos.py`
   already stubs it with a plain fake function (`fake_av_depth`), never touching `ISClipped`.

## Scope

### 1. Rename the method

Rename `ISClipped._av_depth` (`isclipped.py:667-673`, `@lru_cache`-decorated) to
`ISClipped.average_depth`. Update its one internal caller at `isclipped.py:696`
(`report_average`).

### 2. Update the call site

`ijump.py:150` — change `is_processing._av_depth` to `is_processing.average_depth`. No other
change to that call.

### 3. Add direct unit test coverage

`average_depth` has no direct unit test today — only indirect coverage via `report_average`,
itself only exercised through subprocess-level CLI tests. Add one.

`pysamstats.load_coverage` needs a real `pysam.AlignmentFile` — `tests/fake_alignment.py`'s
`FakeAlignment` (only stubs `.references`/`.lengths`) won't work here. Open
`tests/fixtures/tiny.bam` (already used elsewhere in the suite, e.g. `conftest.py`'s
`run_ijump` fixture) via real `pysam.AlignmentFile`, construct an `ISClipped` against it, and
assert `average_depth(chrom, start, stop)` returns the expected mean coverage for a known
region in that fixture.

## Out of scope

- `DetectionResult` / splitting `ISClipped` into an orchestrator + result object (architecture
  review candidate 01).
- Any change to `circos.py`'s signature, `test_circos.py`, or Circos output format.
- Replacing `@lru_cache` with explicit caching.

## Verification

- New unit test for `ISClipped.average_depth` against `tests/fixtures/tiny.bam` passes.
- `pytest` passes from a clean clone.
- `grep -rn "_av_depth"` finds no remaining references outside this ticket's history (i.e. the
  rename is complete — no leftover callers of the old name).

## Done when

- `ISClipped.average_depth` exists; `ISClipped._av_depth` no longer does.
- `ijump.py:150` calls `is_processing.average_depth`.
- A direct unit test for `average_depth` exists and passes.
- `pytest` passes from a clean clone.

## Comments

Implemented 2026-08-13: renamed `ISClipped._av_depth` to `average_depth`
(`isclipped.py:668`), updated its internal caller (`isclipped.py:696`) and the
`ijump.py:150` call site, updated the now-stale name in a `circos.py` comment
(needed to satisfy the ticket's own `grep -rn "_av_depth"` verification step),
and added `tests/test_isclipped.py` exercising `average_depth` against
`tests/fixtures/tiny.bam` via a real `pysam.AlignmentFile`. Full suite passes
except the pre-existing, unrelated `test_read_count_mtx_rejects_invalid_orientation`
failure (references a method that no longer exists, predates this ticket).

**Correction (2026-08-17, review-followups ticket 02):** The `test_read_count_mtx_rejects_invalid_orientation` failure noted above was *not* pre-existing. It was introduced by isclipped-refactor ticket 09's extraction of the read-count-matrix helper to module level in `frequency_estimation.py`, which left no delegating alias on `ISClipped`; on `master` the helper is still an `ISClipped` static method and the test passes there. It is fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.
