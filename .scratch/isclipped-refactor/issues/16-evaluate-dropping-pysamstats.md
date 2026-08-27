# 16 — Drop `pysamstats` from `average_depth`

**What to build:** `ISClipped.average_depth` (`src/ijump/isclipped.py:729`) is
`pysamstats`' only call site in the codebase. Replace it with a pure-`pysam`
implementation and drop the `pysamstats` dependency entirely — or document why
that isn't achievable and keep it.

Two rounds:

- **Round 1 (2026-08-16, closed — do not retry):** swap in pysam's
  `AlignmentFile.count_coverage`. Failed both bars. Full findings in Comments;
  the round-1 body is kept below for the record.
- **Round 2 (2026-08-17, current):** don't use a pileup at all. Port the
  coverage algorithm from `chrcov` — filter reads by SAM flag, sum the lengths
  of reference-consuming CIGAR operations, divide by span length. Scope in
  **Round 2** below.

## Why

`pysamstats` is used exactly once, in `average_depth`. The pysam-only
implementation it replaced is still sitting there, commented out
(`isclipped.py:730-732`):

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

---

## Round 1 — `count_coverage` (closed, rejected)

*Kept verbatim for the record. These steps are **done**; see Comments for the
measured outcome. Do not redo them.*

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
- **Bar:** exact match required.

### Step 2 — Benchmark (only if Step 1 finds an exact match)

- Time both implementations' report-stage cost on a realistically-sized
  dataset with the number of distinct `(chrom, start, stop)` calls a real
  run makes (benchmark through the cached method, not the raw coverage call
  in isolation, since that's what actually runs in production).
- **Bar:** no worse than ~1.5x `pysamstats`' wall-clock time.

**Outcome:** filter semantics reconcilable (`quality_threshold=0`), but exact
match against current production output impossible because of an int32
truncation quirk (see Comments), and `count_coverage` measured 1.85x–2.11x
slower. Both bars failed.

---

## Round 2 — `chrcov`-style CIGAR accumulator

Reference implementation: `chrcov`'s `_accumulate_aligned_bases`
(`~/Documents/WORK/Job_search_2026_04/Jobs/2026-06-15_Harvard/take_home_20260731/harvard_assignment_20260731/src/chrcov/coverage.py`),
an external project of the user's, MIT-licensed and readable but **not** a
dependency — the algorithm gets ported, not imported.

### Why this can work where round 1 didn't

Mean per-base coverage over a span is *sum of aligned bases over span length*.
No per-column pileup is needed to compute it. chrcov exploits that: for each
read, check the SAM flag, then add `sum(length for op, length in
read.cigartuples if op in {0, 7, 8})` (M, `=`, X — the ops that consume both
reference and query, i.e. the only bases sitting on a reference position). It
never touches `SEQ` or base qualities.

That removes both round-1 blockers at the root:

- **Correctness:** there is no base-quality filter to reconcile, because the
  arithmetic never reads base qualities — the same position `pysamstats` is
  already at. Round 1 needed `quality_threshold=0` to get there; here it's free.
- **Performance:** the per-column pileup machinery is exactly what made
  `count_coverage` ~2x slower. A CIGAR sum skips it entirely.

And it's pure `pysam`, which is the actual objective: `pysamstats` leaves the
dependency list, `pytest` runs in a plain venv, tickets 04/05 lose their
conda-only dependency line.

### Correctness bar for round 2 (changed from round 1)

Round 1's bar was byte-for-byte equality with *what `average_depth` returns
today*. That is now known to be unmeetable for a defensible reason: today's
return value is wrong. `average_depth` calls `statistics.mean()` on a
`numpy.int32` array, and `statistics.mean` coerces its result back to the input
element type, truncating every fractional mean to an integer (a region with true
mean 1.5 returns `1`).

**The round-2 bar is exact match against the true, unrounded mean of
`pysamstats`' `reads_all`** — the value today's implementation *should* be
returning. Fixing the truncation is in scope for this ticket (Step 4), and the
`Depth` column will start carrying fractional values as a result. That is a
deliberate, user-visible behavior change, not drift.

### Step 1 — Port the algorithm to a windowed form

chrcov sums whole-read aligned bases in one full-file pass because it reports
whole-contig means. `average_depth(chrom, start, stop)` needs a window, so
reads must be **clipped to `[start, stop)`** — summing whole-read CIGAR lengths
would inflate the mean with bases that overhang the window.

- Use `self.aln.fetch(chrom, start, stop)` and, per read,
  `read.get_blocks()` — pysam already exposes exactly the reference-consuming
  aligned blocks as `(ref_start, ref_end)` pairs — summing
  `max(0, min(block_end, stop) - max(block_start, start))`. Equivalent to
  walking `cigartuples` by hand against a reference cursor, with less to get
  wrong; if `get_blocks()` turns out to merge or split blocks in a way that
  breaks parity, fall back to the explicit CIGAR walk chrcov uses.
- Pick one of two traversal shapes in this step and record why:
  - **(a) per-call `fetch()` over the window.** Simplest, keeps `lru_cache`
    and the current call shape, cost scales with reads overlapping the window.
    **Start here.**
  - **(b) one full-file pass building a per-contig `+1/-1` diff array
    (numpy) → `cumsum` → prefix sums**, making every subsequent
    `average_depth` an O(1) lookup. Costs ~4-8 bytes per reference base —
    fine for the few-Mb bacterial genomes iJump targets, unbounded on larger
    references. Only reach for this if (a) misses the Step 3 bar.

### Step 2 — Reconcile semantics against `pysamstats`

Round 1 died here, on a difference nobody predicted. The divergences are
different this time, and each must be *tested*, not assumed:

- **Denominator.** `pysamstats.load_coverage(..., truncate=True)` is called
  with default `pad=False`, which emits a row only for reference positions the
  pileup actually produced — so today's `mean(c.reads_all)` averages over
  *covered* positions, not over `stop - start`. A CIGAR accumulator naturally
  divides by window length. Test a window containing zero-coverage positions,
  determine which denominator each uses, and match `pysamstats`' behavior
  exactly (or, if `pad=False` turns out not to skip positions here, document
  that and use window length).
- **Supplementary reads.** htslib's pileup default filter excludes
  UNMAP|SECONDARY|QCFAIL|DUP but **not** SUPPLEMENTARY; chrcov's default mask
  (3844) *does* exclude supplementary. Porting the mask verbatim would silently
  drop reads `pysamstats` counts. Keep supplementary reads (chrcov's
  `--include-supplementary` behavior) unless a test shows otherwise. This is
  not cosmetic: iJump is a transposon-insertion caller, and supplementary
  alignments of clipped reads are precisely the reads it exists to find.
- **Deletions and reference skips.** A read spanning a `D` or `N` still
  occupies those pileup columns (`is_del`/`is_refskip`), and `reads_all` counts
  every read in a column. A sum over M/`=`/X alone does not count them. Build a
  read with an internal deletion inside the window, measure whether `reads_all`
  includes it, and either add `D` (op 2) and `N` (op 3) to the counted ops or
  document the difference — if it can't be closed, that's a round-2 stop the
  same way round 1 stopped.
- **Base quality.** Nothing to do — neither implementation filters on it.
- **`max_depth=300000`** on the current call caps pileup depth; the accumulator
  has no cap. Divergence only above 300000x coverage. Note it, don't engineer
  for it.
- **Test fixtures.** Extend
  `tests/test_average_depth_pysamstats_vs_count_coverage.py` (rename it, or add
  a sibling — it's now covering more than count_coverage): its synthetic-BAM
  harness and 80x/200kb generator are reusable, but its fixtures have no
  supplementary-flagged reads, no reads with internal deletions, and no
  zero-coverage windows. Add all three.
- **Bar:** 0 mismatches against the true float mean of `reads_all` across the
  random-window sweep, at the same scale round 1 used (50-region and 300-region
  sweeps).

### Step 3 — Benchmark

- Reuse round 1's benchmark harness verbatim: synthetic 80x/200kb BAM, 300
  distinct random windows of typical IS-flank sizes, timed through the cached
  `average_depth`. Round 1's baseline numbers (pysamstats ~1.42–1.50s) are
  directly comparable.
- **Bar:** no worse than ~1.5x `pysamstats`. The expectation is *faster* than
  pysamstats, since there's no per-column pileup on either side of the
  comparison — so a ratio near or above 1.5x is a signal that traversal shape
  (b) should be tried before this is called a failure, not an immediate stop.

### Step 4 — Swap, and fix the truncation

- Replace `average_depth`'s body with the validated accumulator, returning a
  `float`. Drop the `pysamstats` import, the now-unused `statistics.mean`
  import, and the commented-out `count_coverage` lines.
- **Audit downstream for integer-depth assumptions.** The `Depth` column now
  carries fractional values where it previously truncated. Check
  `region_summary.py`, `frequency_estimation.py`, and the report writers for
  formatting, dtype, or division that assumes an int. Frequency estimation
  divides by depth, and truncating to a *smaller* integer inflated frequencies
  — reported frequencies will shift slightly downward. Record this in Comments
  and flag it for the README/changelog as a user-visible change.
- Remove `pysamstats` wherever it's declared at implementation time:
  `README.md`'s conda-install section, `environment.yml`, `meta.yaml`,
  `pyproject.toml`, `Dockerfile` (tickets 04/05 own those files' content —
  coordinate rather than rewrite them wholesale).
- Update `src/ijump/circos.py:14`'s comment if it still names pysamstats.

## Out of scope

- Any other coverage/depth calculation elsewhere in the codebase — grepped
  at grilling time, `pysamstats` has exactly one call site, so this is a
  single-method change, not a codebase-wide sweep.
- Changing `average_depth`'s signature, caching behavior, or callers (the
  return *type* changing from truncated-int to float is intended and is not
  a signature change).
- Porting anything else from chrcov: its read-accounting/provenance log, its
  CLI, its CRAM reference handling, its whole-contig reporting. The algorithm
  is the deliverable; chrcov does not become a dependency.
- Adding a MAPQ filter. chrcov has `--min-mapq`; `pysamstats`' pileup applies
  none, so introducing one would break parity by construction. Leave it out.
- Fixing the `lru_cache`-on-instance-method leak flagged at
  `isclipped.py:724` (ticket 06 already owns that note).

## Verification

- Parity test committed regardless of outcome, covering supplementary reads,
  internal deletions, zero-coverage windows, and a random-window sweep.
- If swapped: `pytest` passes **in a plain venv with `pysamstats` not
  installed** — this is the actual prize, and the only check that proves the
  dependency is really gone. Grep confirms no `pysamstats` import anywhere in
  `src/ijump`.
- If swapped: manual real-sample run (reused from prior tickets) compared
  against a pre-swap run. `Depth` values must differ from the old ones *only*
  by the intended de-truncation — assert the relationship between old and new
  values in the comparison rather than assuming which rounding rule
  `statistics.mean`'s coercion applied.
- Benchmark ratio documented in Comments regardless of outcome.
- If not swapped: Comments document the measured gap that failed the bar, and
  `pysamstats` stays in tickets 04/05's dependency lists.

**Blocked by:** None — can start immediately. Note for whoever picks up
packaging tickets 04/05: their `pysamstats` dependency line is contingent
on this ticket's outcome — check its `Status`/`Comments` before finalizing
`environment.yml`/Dockerfile dependency lists.

**Status:** done

Round 1 (closed, do not redo):

- [x] Characterization test comparing `count_coverage` vs `pysamstats.load_coverage` committed.
- [x] Exact-match bar evaluated and documented — failed.
- [x] Performance measured anyway (1.85x–2.11x) and documented.

Round 2:

- [x] Windowed CIGAR accumulator implemented, traversal shape chosen and justified.
- [x] Semantics reconciled and tested: denominator (`pad=False` covered-positions vs window length), supplementary reads, deletions/ref-skips.
- [x] Parity bar evaluated against the true float mean of `reads_all` and documented.
- [x] Benchmark run against round 1's harness; ratio documented.
- [x] If both bars met: `average_depth` swapped, truncation fixed, downstream integer-depth assumptions audited.
- [x] If swapped: `pysamstats` import and dependency declarations removed; `circos.py:14` comment updated.
- [x] `pytest` passes in a plain venv with `pysamstats` uninstalled.
- [x] Frequency-estimation impact of the de-truncation documented for the changelog/README.
- [x] Outcome (swapped or not, with numbers) documented in Comments for tickets 04/05.

## Comments

### Reopened 2026-08-17 — round 2, after grilling on a chrcov-style algorithm

Round 1 closed this on the premise that the only pure-pysam alternative was
`count_coverage`. It isn't. `chrcov` computes mean coverage with no pileup at
all — flag filter, then sum reference-consuming CIGAR op lengths, divide by
span. Grilling established that this attacks both round-1 failures at the root
rather than working around them: no base-quality filter exists to reconcile
(the arithmetic never reads `SEQ` or quality strings), and the per-column
pileup that made `count_coverage` ~2x slower is skipped entirely.

Two decisions from that session, both changing the ticket's shape:

1. **The correctness bar moves from "byte-for-byte with today's output" to
   "exact match against the true float mean."** Round 1 correctly identified
   that current output truncates fractional means to integers (the
   `statistics.mean` + `numpy.int32` coercion quirk), and correctly treated
   that as a stop under its own bar — but the conclusion should have been that
   the bar was wrong, not that the swap was. Fixing the truncation is now in
   scope here rather than deferred, since any correct replacement fixes it
   incidentally.
2. **New parity risks replace the old ones**, and they're the real work of
   round 2: the `pad=False` denominator question (today's mean may be over
   covered positions, not window length), supplementary-read handling (chrcov's
   default mask excludes them, htslib's pileup default does not — and clipped
   supplementary alignments are iJump's core signal), and deletion/ref-skip
   columns (counted by `reads_all`, not by a sum over M/`=`/X). Each is a
   candidate round-2 stop and each must be measured, not reasoned about.

Round 1's test harness and benchmark rig are reusable as-is; its fixtures are
not sufficient (no supplementary flags, no internal deletions, no zero-coverage
windows).

### Round 1 — evaluated 2026-08-16. Decision: do not swap to `count_coverage`.

#### Step 1 — correctness: filter parity achievable, but not exact-match overall

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

*(Round-2 note: this reasoning is why the bar itself was revisited — see the
reopening comment above. The measurement stands; the conclusion drawn from it
no longer does.)*

#### Step 2 — performance: measured anyway, also fails the bar

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

#### Round 1 net

Both bars fail independently. `pysamstats` is not dropped *via `count_coverage`*.
The characterization test
(`tests/test_average_depth_pysamstats_vs_count_coverage.py`) is committed
regardless, per the ticket's own requirement, and documents both findings —
including the numpy-int32-truncation quirk in current production behavior.

Full `pytest` run: 45 passed, 1 pre-existing failure unrelated to this ticket
(`test_read_count_mtx_rejects_invalid_orientation` — references a method,
`ISClipped._read_count_mtx`, that no longer exists; predates this ticket, noted
in ticket 12's comments too).

**Correction (2026-08-17, review-followups ticket 02):** The `test_read_count_mtx_rejects_invalid_orientation` failure noted above was *not* pre-existing. It was introduced by isclipped-refactor ticket 09's extraction of the read-count-matrix helper to module level in `frequency_estimation.py`, which left no delegating alias on `ISClipped`; on `master` the helper is still an `ISClipped` static method and the test passes there. It is fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.

### Round 2 — evaluated and shipped 2026-08-17. Decision: swap to the CIGAR/span accumulator.

#### Step 1 — traversal shape: (a), per-call `fetch()`

Went with shape (a) from the ticket (per-call `fetch()` over the window, kept
under `lru_cache` with the current call shape) since it cleared the Step 3
benchmark bar comfortably; shape (b)'s prefix-sum array was never needed.

The implementation (`ISClipped.average_depth`, `isclipped.py`) does not use
`get_blocks()` as the ticket's Step 1 first suggested — verified directly
that `get_blocks()` excludes deletion/ref-skip spans (a read `5M 5D 5M` at
pos 10 returns blocks `[(10, 15), (20, 25)]`, dropping `[15, 20)`), which
would have undercounted relative to `pysamstats`. Simpler than either
`get_blocks()` or an explicit CIGAR walk: a read's
`[reference_start, reference_end)` span is *already* contiguous coverage for
`pysamstats`' purposes, because every reference-consuming CIGAR op (M/=/X,
but also D/N) sits inside that span and only S/I/H (which don't consume the
reference) fall outside it. So the accumulator just clips
`(reference_start, reference_end)` to `(start, stop)` per read and sums the
overlap into a small `numpy` array — no CIGAR walking, no `get_blocks()`.

#### Step 2 — semantics reconciled, all three confirmed by direct measurement (not assumed)

Using a live `pysam`/`pysamstats` install recovered from a conda-env backup
(`~/miniconda3_envs_backup/ijump-verify`, pysam 0.23.0 / pysamstats 1.1.2 —
neither is installed in this repo's own `.venv` or the base conda env, so
this was necessary groundwork before any of Step 2 could be measured):

- **Denominator.** Confirmed `pad=False` (the default, and what production
  used) emits a `pysamstats` row *only* for positions some read actually
  covers — a window `[0, 30)` with one read spanning `[10, 20)` returns 10
  rows (`pos` 10-19), not 30. `average_depth`'s denominator is therefore
  "count of covered positions in the window", not window length — see
  `test_denominator_is_covered_positions_not_window_length`.
- **Supplementary reads.** Confirmed by direct test: a region with one
  normal read and one supplementary-flagged (0x800) read at the same
  position returns `reads_all` mean 2.0, i.e. htslib's pileup default (which
  `pysamstats` sits on) does *not* exclude supplementary reads. Kept them
  counted, matching the ticket's expectation — see
  `test_supplementary_reads_are_counted`.
- **Deletions/ref-skips.** Confirmed by direct test: a read `5M 5D 5M` (or
  `5M 5N 5M`) at pos 10 produces uniform depth-1 `reads_all` across all 15
  positions `[10, 25)`, including the 5 positions inside the D/N span. Since
  the accumulator uses the read's full `reference_start:reference_end` span
  (see Step 1), this falls out for free without special-casing D/N ops — see
  `test_internal_deletion_counts_as_covered`, `test_ref_skip_counts_as_covered`.
- **Base quality.** As anticipated, nothing to reconcile — neither side
  filters on it.
- **`max_depth=300000`.** Not engineered for, as scoped.

`tests/test_average_depth_pysamstats_parity.py` (renamed from
`test_average_depth_pysamstats_vs_count_coverage.py`) carries all of this:
round-1's tests kept (marked `test_round1_*`, historical record, still
true), a fixed version of the truncation test (now asserts production
returns the untruncated float), and new round-2 tests for each semantic
axis above plus a 300-window sweep (`n_reads` ≈ 80x/200kb, with
supplementary/duplicate/secondary flags and internal deletions/ref-skips
mixed into ~12% of reads) comparing the real, unmodified
`ISClipped.average_depth` against `pysamstats`' true float mean.

**Result: 0/300 mismatches.** Every `pysamstats`-dependent test in the file
calls `pytest.importorskip("pysamstats")` itself (not at module level), so
the file collects and passes (with those tests skipped, not failed) in a
plain venv with `pysamstats` absent, and runs the full parity comparison in
an environment that has it.

#### Step 3 — benchmark: faster than pysamstats, not just within bar

Reused round 1's harness shape (synthetic 80x/200kb BAM, 300 distinct random
windows of realistic IS-flank sizes, timed through the cached
`average_depth`). Three runs, same script structure as round 1's:

| run | pysamstats | accumulator | ratio |
|-----|-----------|--------------|-------|
| 1   | 0.718s    | 0.465s       | 0.65x |
| 2   | 0.798s    | 0.466s       | 0.58x |
| 3   | 0.721s    | 0.462s       | 0.64x |

**~1.5-1.7x faster than `pysamstats`**, not just under the ~1.5x-slower bar —
matching the ticket's prediction that skipping the per-column pileup
entirely (not just tuning it) would flip the direction of round 1's
regression.

#### Step 4 — swapped, truncation fixed, downstream audited

- `ISClipped.average_depth` (`isclipped.py`) now returns a `float`
  unconditionally (raises `statistics.StatisticsError` for a window with
  zero covered positions, matching what `mean()` on an empty `pysamstats`
  array already did in production). `pysamstats` import,
  `from statistics import mean`, and the commented-out `count_coverage`
  lines are gone; replaced with `from statistics import StatisticsError`
  (used for the empty-window case) and a module-level `_COVERAGE_EXCLUDE_FLAGS`
  constant documenting the flag mask.
- **Downstream audit:** `region_summary.py`'s `report_average` (the only
  caller that feeds `average_depth`'s return value into the `Depth` column
  in `--estimation_mode average`) does plain float division/`round(...)` —
  no int assumption. `circos.py`'s two `average_depth_fn` call sites do
  `depth > 0` comparisons, float division, and `str(depth)` — all
  float-safe. `frequency_estimation.py`'s own `Depth` column (precise mode)
  is unrelated — it sums `N_*` read counts, never calls `average_depth`.
  `tests/test_isclipped.py`'s existing `== 2` assertion against
  `tests/fixtures/tiny.bam` still passes (`2 == 2.0`).
- **Manual real-sample comparison:** `Test/TCL_A1/Slice.bam` (real data, not
  synthetic), 200 windows anchored on real read start positions, old
  (`pysamstats`+`statistics.mean`) vs. new (production `average_depth`)
  compared directly. **0 relationship violations**: for every non-empty
  window, `int(new_value) == old_value` (new is old's untruncated float —
  e.g. old `415` → new `415.84293193717275`); for the (this slice had none
  in the anchored sample, but verified separately) zero-coverage case, new
  raises where old's `mean()` would have raised on an empty array.
- **Dependency removal:** `pyproject.toml`, `environment.yml`, `meta.yaml`,
  `README.md`'s dependency list/conda section, and the `Dockerfile`'s
  `pysamstats` install step all had their `pysamstats` declarations
  removed. `circos.py:14`'s comment (`"needs a live BAM handle via
  pysamstats"`) updated to say `pysam`. Left as-is (out of this ticket's
  scope, flagged instead): the `Dockerfile`'s Python 3.8 pin / `pysam==0.15.4`
  / `build-essential` / `zlib1g-dev` machinery, which existed solely to
  satisfy `pysamstats`' `pysam<0.16` pin and could now be relaxed to a
  modern `pysam` wheel on a current Python base — a job for whoever next
  touches packaging ticket 05 (Docker), not re-architected here.
- **Plain-venv verification (the actual prize):** built a clean `uv venv
  --python 3.11` at `/tmp/ijump-plain-venv` (this repo's own `.venv` is
  Python 3.13, which hits an unrelated `mypy==2.3.1` `requires-python>=3.10`
  resolution conflict via `uv sync`'s dependency groups — sidestepped by
  installing runtime deps directly rather than via `uv sync`), installed
  `pysam` + `biopython==1.79` + `pandas<3` + `numpy` + `scikit-learn` +
  `pytest` + `ijump` (editable, `--no-deps`) with plain `uv pip install`,
  confirmed `pysamstats` is not importable there, and ran the full suite:
  **50 passed, 3 skipped** (the 3 are `pysamstats`-gated parity tests,
  correctly skipped rather than failed). `grep -rn pysamstats src/ijump`
  finds only comments (no imports) in the changed files.
- **Frequency-estimation impact documented:** added an "Unreleased" entry
  at the top of `README.md` (ahead of the existing v1.0.4 changelog note)
  explaining that `--estimation_mode average`'s `Depth` column now carries
  fractional values, and `Frequency` (which divides by `Depth`) will shift
  slightly *downward* relative to old runs, since the old truncated
  (smaller) `Depth` inflated it. Also updated the nearby "tests additionally
  require `pysam` and `pysamstats`" sentence in the dev-setup section; left
  the larger historical `uv`/Docker sections documenting the now-resolved
  `pysam<0.16` conflict as-is (still accurate as history, flagged with a
  pointer note) rather than rewriting them wholesale.

#### Post-implementation review (code-review skill, standards + spec axes)

Both axes ran clean. Standards found no hard violations (a `Shotgun Surgery`
judgement call for the dependency removal touching 8 files was noted but is
the inherent shape of this kind of ticket, and a mild `Primitive Obsession`
call on the raw-int `_COVERAGE_EXCLUDE_FLAGS` bitmask, both low severity).
Spec found one real gap, since fixed: Step 2's bar calls for parity sweeps
"at the same scale round 1 used (50-region and 300-region sweeps)", but the
round-2 sweep against the actual swapped `average_depth` only covered 300
windows, not 50. Fixed by parametrizing
`test_production_matches_pysamstats_across_realistic_windows_with_gaps_and_supplementary`
over `n_windows in (50, 300)` — both scales now pass with 0 mismatches
against the real production method. Spec also flagged the
`Dockerfile`/`README.md` edits as going slightly past literal
"remove pysamstats declarations" (added advisory prose about Dockerfile's
now-unnecessary Python 3.8 pin) — judged in-spirit with the ticket's own
"coordinate rather than rewrite them wholesale" instruction and left as-is,
but noted here for whoever picks up ticket 05.

#### Round 2 net

Both bars pass, independently and by a comfortable margin (0/300 parity
mismatches; 1.5-1.7x faster, not just under the 1.5x-slower ceiling).
`pysamstats` is fully dropped from the runtime dependency surface;
`environment.yml`/`meta.yaml`/`Dockerfile` (tickets 04/05) no longer need to
carry it. `lint`/`mypy`/full `pytest` all pass, including in a plain venv
with `pysamstats` uninstalled — the condition tickets 04/05 were blocked on
checking before finalizing their own dependency lists.
