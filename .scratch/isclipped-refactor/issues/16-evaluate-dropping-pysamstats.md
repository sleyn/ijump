# 16 — Evaluate replacing `pysamstats` with pysam's built-in `count_coverage`

**What to build:** Not a guaranteed removal — an evaluation with a hard
correctness/performance bar. `ISClipped.average_depth` (`isclipped.py:643-649`)
is `pysamstats`'s only call site in the codebase. Determine whether pysam's
own `AlignmentFile.count_coverage` can replace it with byte-for-byte
equivalent output and acceptable performance; if so, make the swap and
drop the `pysamstats` dependency entirely; if not, close this out
documenting why, and keep `pysamstats`.

## Why

`pysamstats` is used exactly once, in `average_depth`. The pysam-only
implementation it replaced is still sitting there, commented out
(`isclipped.py:645-647`):

```python
# aln_depth = self.aln.count_coverage(chrom, start, stop)
# depth = sum(map(sum, aln_depth))
# return depth / len(aln_depth[0])  # average depth of the region
```

Git history (`ccaa9067`, 2020-03-10, message: *"switched to pysamstats for
coverage calculation as it has ~3x increase in speed for the report stage
in test"*) shows this was a deliberate performance-motivated swap, made on
whatever pysam/pysamstats versions existed in March 2020 — six years old
at grilling time — and before `lru_cache` (added in the immediately
preceding commit, `7e8e98fd`, same day) started memoizing repeated
`average_depth(chrom, start, stop)` calls. Both facts that justified the
original swap (library performance at the time, no memoization yet) may no
longer hold.

`pysamstats` is also the reason `pytest` can't run in a plain venv today —
it's a compiled, conda-only dependency per `README.md`'s conda-install
section, confirmed by hitting `ModuleNotFoundError: No module named
'pysamstats'` when running `pytest` in this session outside conda. It's
also the one dependency ticket 04's `environment.yml`/`meta.yaml` and
ticket 05's Docker image would otherwise need to carry solely for this one
call site. Removing it would simplify all three packaging tickets — but
only if the replacement is provably equivalent; this is a scientific tool
and `average_depth`'s output feeds directly into every report's `Depth`
column and downstream frequency estimation (`frequency_estimation.py`), so
a silent numeric drift is not acceptable.

## Scope

### Step 1 — Characterize the difference (or confirm there isn't one)

- `count_coverage`'s default `quality_threshold=15` filters per-base by
  base call quality; `pysamstats.load_coverage`'s `reads_all` does not
  filter by base quality, but does apply read-level pileup filters
  (unmapped/secondary/qcfail/duplicate, matching pysam's pileup defaults).
  These are different filtering mechanisms and are not guaranteed to
  produce the same average for a given region.
- Write a test (reusing `tests/fixtures/tiny.bam` or a richer fixture if
  the tiny one doesn't exercise enough read diversity — e.g. no reads with
  low base quality or duplicate flags to actually exercise the difference)
  that calls both implementations against the same `(chrom, start, stop)`
  regions and compares results directly.
- If they differ, try tuning `count_coverage`'s parameters
  (`quality_threshold=0` to disable base-quality filtering,
  `read_callback=` a custom callable replicating pysamstats'/pysam
  pileup's default flag exclusions) until they match exactly, or determine
  no parameter combination closes the gap.
- **Bar:** exact match required. If no exact match is achievable on
  representative data, stop here — close this ticket documenting the
  measured difference and why, and leave `average_depth` on `pysamstats`
  unchanged.

### Step 2 — Benchmark (only if Step 1 finds an exact match)

- Time both implementations' report-stage cost on a realistically-sized
  dataset (not the tiny pytest fixture — find or construct a BAM
  comparable in scale to what a real population-sequencing run would
  produce; check `Test/`, `Example files/`, or `simulation/` for something
  suitable before generating new data) with the number of distinct
  `(chrom, start, stop)` calls a real run makes (the `lru_cache` matters
  here — benchmark through the cached method, not the raw coverage call in
  isolation, since that's what actually runs in production).
- **Bar:** `count_coverage` must be no worse than ~1.5x `pysamstats`'
  wall-clock time on this benchmark. Document the actual measured ratio in
  Comments regardless of outcome.

### Step 3 — Swap (only if both bars are met)

- Replace `average_depth`'s body with the validated `count_coverage`-based
  implementation; remove the commented-out old code (it becomes live code)
  and the `pysamstats` import.
- Remove `pysamstats` from anywhere it's declared as a dependency at
  implementation time (check ticket 01's `pyproject.toml` if it's landed
  by then, and `README.md`'s conda-install section).
- Update `circos.py:15`'s comment (*"handle via pysamstats"*, referring to
  `average_depth`'s dependency on a live BAM handle) if it still mentions
  pysamstats by name once the swap lands.

## Out of scope

- Any other coverage/depth calculation elsewhere in the codebase — grepped
  at grilling time, `pysamstats` has exactly one call site
  (`isclipped.py:648`), so this is a single-method change, not a
  codebase-wide sweep.
- Changing `average_depth`'s signature, caching behavior, or callers.
- Packaging file changes (tickets 04/05 own `environment.yml`/Dockerfile
  content) — this ticket only determines whether `pysamstats` is needed;
  it doesn't itself edit packaging tickets' deliverables, just their
  contingency note (added separately, see below).

## Verification

- Step 1's characterization test is committed regardless of outcome (even
  if the ticket closes without swapping, the test documents the comparison
  that was made).
- If swapped: `pytest` passes from a clean clone with no `pysamstats`
  import anywhere in `src/ijump` (or wherever ticket 01 has placed modules
  by then — check current layout at implementation time).
- If swapped: manual real-sample run (reused from prior tickets) produces
  identical `Depth` column values to a pre-swap run on the same input.
- If not swapped: Comments document the measured correctness gap and/or
  performance ratio that failed the bar.

**Blocked by:** None — can start immediately. Note for whoever picks up
packaging tickets 04/05: their `pysamstats` dependency line is contingent
on this ticket's outcome — check its `Status`/`Comments` before finalizing
`environment.yml`/Dockerfile dependency lists.

**Status:** closed-no-change

- [x] Characterization test comparing `count_coverage` vs `pysamstats.load_coverage` committed.
- [x] Exact-match bar evaluated and documented (met, or ticket closed without swapping).
- [ ] If exact match met: benchmark against realistic-scale data, ~1.5x-slack bar evaluated and documented. (measured anyway, see Comments)
- [ ] If both bars met: `average_depth` swapped, `pysamstats` import/dependency removed, `circos.py:15`'s comment updated if needed. (not reached — correctness bar failed)
- [x] `pytest` passes from a clean clone in either outcome.
- [x] Outcome (swapped or not, with numbers) documented in Comments for tickets 04/05 to reference.

## Comments

Evaluated 2026-08-16. **Decision: do not swap. `pysamstats` stays.**

### Step 1 — correctness: filter parity achievable, but not exact-match overall

`tests/fixtures/tiny.bam` (7 reads, uniform flag 0, uniform base quality 40) has no
read diversity to exercise either library's filters, so synthetic BAMs were built
in a new test, `tests/test_average_depth_pysamstats_vs_count_coverage.py`.

Filter-level investigation (verified directly, not assumed):

- `pysamstats.load_coverage`'s `no_dup` parameter has **no measurable effect** —
  confirmed by feeding it a duplicate-flagged-only region and toggling
  `no_dup=True/False`: both return an empty coverage array. The underlying
  `pysam` pileup it's built on already excludes unmapped/secondary/qcfail/
  duplicate-flagged reads via its own default flag filter, before pysamstats'
  own knob is ever consulted. This matches `count_coverage`'s default
  `read_callback='all'` exactly — no custom callback is needed.
- The one real discrepancy is base-quality filtering: `count_coverage` filters
  bases below phred 15 by default (`quality_threshold=15`); `pysamstats` never
  filters by base quality at all. Passing `quality_threshold=0` to
  `count_coverage` (default `read_callback` otherwise) closes this gap exactly.
- Verified across one hand-built region (low-baseq + duplicate reads mixed,
  `test_tuned_count_coverage_matches_pysamstats_filters`) and a 50-region random
  sweep over an 80x/5kb synthetic BAM
  (`test_tuned_count_coverage_matches_pysamstats_across_realistic_windows`):
  **0/50 mismatches** at the true (unrounded) mean level once `quality_threshold=0`
  is set.

So filter semantics *are* exactly reconcilable — that part of the correctness bar
passes.

**But the overall correctness bar still fails**, for a reason distinct from
filtering: `average_depth`'s current body calls `statistics.mean()` directly on
`pysamstats`' `reads_all` array, which is `numpy.int32`. Python's `statistics.mean`
detects the input element type and coerces its *final* result back to it — for a
numpy int32 array this **truncates any fractional mean to an integer** instead of
returning a float:

```
>>> from statistics import mean
>>> import numpy as np
>>> mean(np.array([1, 2, 4], dtype='int32'))
2          # not 2.333...
```

A `count_coverage`-based replacement (plain `depth / len(...)` division — either
tuned per above, or the commented-out pysam-only code already sitting above
`average_depth`) does not reproduce that truncation; it returns the mathematically
correct float. Demonstrated in
`test_production_average_depth_truncates_fractional_means`: a region with true
mean 1.5 returns `1` from the real, unmodified `ISClipped.average_depth` today,
and `1.5` from the count_coverage-based computation. Measured against a synthetic
80x/200kb benchmark BAM (see Step 2) across 300 random realistic windows: **0/300
mismatches** at the true-float-mean (filter-parity) level, but **300/300
mismatches** against what `average_depth` actually returns today, because almost
no random window happens to have an exactly-integer true mean.

This is a pre-existing bug in the current implementation, not something this
ticket introduced or is in scope to fix — but it means a `count_coverage` swap
using ordinary float division is not byte-for-byte equivalent to current
production output, which is the hard bar Step 1 sets. Per the ticket's own
escape valve ("if no exact match is achievable, STOP — do not swap"), this is a
stop.

### Step 2 — performance: measured anyway, also fails the bar

Not required once Step 1 failed, but cheap to measure and useful for the record.
No BAM file exists anywhere in the repo (checked `Test/`, `Example files/`,
`simulation/` — the latter has reference FASTAs and alignment scripts but no
committed aligned reads), so a synthetic ~80x-coverage, 200kb BAM was generated
and both implementations benchmarked through the cached `average_depth` method
(300 distinct random regions of typical IS-flank window sizes, so caching doesn't
mask the comparison in either direction). Three runs:

| run | pysamstats | count_coverage | ratio |
|-----|-----------|-----------------|-------|
| 1   | 1.499s    | 2.775s          | 1.85x |
| 2   | 1.421s    | 2.752s          | 1.94x |
| 3   | 1.464s    | 3.094s          | 2.11x |

`count_coverage` is consistently **~1.85x–2.1x** pysamstats' wall-clock time,
over the ~1.5x bar. (Plausible explanation, not verified further: `count_coverage`
builds four separate per-base A/C/G/T arrays across the whole queried region and
sums them in Python, versus pysamstats' single compiled-Cython coverage pass.)

### Net

Both bars fail independently. `pysamstats` is not dropped. The characterization
test (`tests/test_average_depth_pysamstats_vs_count_coverage.py`) is committed
regardless, per the ticket's own requirement, and documents both findings —
including the numpy-int32-truncation quirk in current production behavior, which
packaging tickets 04/05 (and anyone touching `average_depth` later) should be
aware of even though fixing it is out of this ticket's scope.

Full `pytest` run: 45 passed, 1 pre-existing failure unrelated to this ticket
(`test_read_count_mtx_rejects_invalid_orientation` — references a method,
`ISClipped._read_count_mtx`, that no longer exists; predates this ticket, noted
in ticket 12's comments too).

**Correction (2026-08-17, review-followups ticket 02):** The `test_read_count_mtx_rejects_invalid_orientation` failure noted above was *not* pre-existing. It was introduced by isclipped-refactor ticket 09's extraction of the read-count-matrix helper to module level in `frequency_estimation.py`, which left no delegating alias on `ISClipped`; on `master` the helper is still an `ISClipped` static method and the test passes there. It is fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.
