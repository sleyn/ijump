# 16 — Evaluate dropping `pysamstats` from `ISClipped.average_depth`

Status: closed-no-change
Blocked by: —

## Why

`ISClipped.average_depth` (`isclipped.py:643-649`) computes average depth over a
region via `pysamstats.load_coverage`. Directly above it sits commented-out code
using pysam's own `self.aln.count_coverage(chrom, start, stop)` instead — an older
implementation that predates the `pysamstats` dependency. `pysamstats` is a real
packaging cost (it pulls in its own compiled extension, pins older pysam/numpy
versions transitively, and is the kind of dependency the packaging tickets (04/05)
would rather not carry if the tool doesn't actually need it). This ticket asks
whether the pysam-only `count_coverage` implementation can exactly replace it.

Two hard bars, in order:

1. **Correctness (exact match required).** `count_coverage`-based average depth
   must exactly match `pysamstats.load_coverage`-based average depth over the same
   `(chrom, start, stop)` regions — including on regions with read diversity
   `tests/fixtures/tiny.bam` doesn't have (low base quality, duplicate flag), since
   `count_coverage`'s default `quality_threshold=15` base-quality filter and
   pysamstats' read-level pileup filters don't obviously agree. Tuning
   `count_coverage`'s params is fair game. If no exact match is achievable: stop,
   don't swap, document the gap, close the ticket.
2. **Performance (only if #1 passes).** Benchmarked through the cached
   `average_depth` method (not the raw coverage call) on realistically-sized data,
   `count_coverage` must be no worse than ~1.5x pysamstats' wall-clock time.
3. **Swap (only if both pass).** Replace `average_depth`'s body, remove the
   commented-out old code and the `pysamstats` import, remove `pysamstats` from
   README.md's conda-install section, update `circos.py:15`'s comment if it names
   pysamstats.

## Scope

- Write a characterization test comparing the two implementations on regions
  exercising both libraries' filtering mechanisms, using `tests/fixtures/tiny.bam`
  or a richer fixture built for this purpose if `tiny.bam` isn't diverse enough.
- Benchmark both implementations through `average_depth` on realistically-sized
  data (check `Test/`, `Example files/`, `simulation/` before generating new data).
- Swap only if both bars are cleared; otherwise stop and document why.

## Out of scope

- Any other coverage/depth call site (there is only this one).
- `average_depth`'s signature, caching (`@lru_cache`), or callers.
- `environment.yml` / `Dockerfile` — owned by the packaging tickets (04/05).

## Verification

- The characterization test is committed either way.
- If swapped: `pytest` passes, `grep -rn pysamstats` finds nothing in source.
- If not swapped: the measured correctness gap and/or performance ratio that
  failed the bar is documented below, for packaging tickets 04/05 to reference.

## Done when

- [x] Characterization test comparing `count_coverage` vs `pysamstats.load_coverage`
      average depth exists and is committed, using a fixture with read diversity
      (low base quality, duplicate flag) `tiny.bam` lacks.
- [x] Correctness bar evaluated: exact match or documented gap.
- [x] Performance bar evaluated (only meaningful if correctness passed, but
      measured regardless here since it's cheap and worth recording).
- [ ] ~~`average_depth` swapped to `count_coverage`, old code/import removed,
      README.md and `circos.py:15` updated~~ (not reached — correctness bar failed).
- [x] `pytest` passes (pre-existing unrelated failure aside, same as before this
      ticket).

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
detects the input element type and coerces its *final* result back to it —for a
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
