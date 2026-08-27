# 01 — Handle zero-coverage regions in average_depth without crashing

Status: ready-for-agent
Blocked by: —

## The bug

`ISClipped.average_depth` (src/ijump/isclipped.py:735) raises an uncaught
`statistics.StatisticsError` when the queried region `[start, stop)` has zero
covered positions (e.g. a fully-deleted annotated region). Nothing catches
`StatisticsError` anywhere in `isclipped.py`, `region_summary.py`, or
`circos.py`, so this crashes the entire run instead of degrading gracefully
for that one region.

`region_summary.report_average` only calls `average_depth` for rows where
`count > 0` (a junction was actually found there), so this is reachable in
practice: a region can have junction-supporting reads and still have zero
total coverage in `[start, stop)`, e.g. when the annotated window is small
relative to the deletion, or when the junction-supporting reads' aligned
spans don't happen to overlap the window at all.

## Fix

1. `average_depth` returns `0.0` instead of raising when zero positions in
   `[start, stop)` are covered.
2. `region_summary.report_average`'s `Frequency` column (region_summary.py:85-92)
   currently divides by `Depth * (1 - min_match / av_read_len)`. With
   `Depth = 0` this becomes `inf` via ordinary float division rather than
   crashing. Clip `Frequency` to `NaN` when `Depth == 0` instead of leaving it
   as `inf`.

## Why these are the right values (context from grilling session)

- `average_depth`'s covered-only averaging (excluding zero-coverage positions
  from the mean, matching pysamstats' `pad=False` semantics) is itself
  correct and out of scope for this ticket. Commit `d671080` (2021) tried the
  alternative — whole-region averaging with `pad=True`, including
  zero-coverage positions in the mean — and abandoned it ("Seems not working
  well") after hitting real deletion regions; the author confirmed this was
  the reason during the grilling session that produced this ticket. That
  commit is not an ancestor of current HEAD.
- `0.0` (not raising, not `NaN`) was chosen for `average_depth`'s return value
  because `circos.py:138` already has a dormant `if depth > 0:` guard around
  its histogram cutoff logic that is currently dead code (the crash happens
  before it's ever reached) — this independently confirms `0.0` aligns with
  pre-existing intent elsewhere in the codebase.
- `NaN` (not `inf`, not `1.0`) was chosen for the resulting `Frequency`
  because a fully-deleted region has no background reads to normalize
  against, so the insertion's population frequency genuinely cannot be
  computed there — `NaN` reads as "not computable," which is more honest
  than asserting a bounded frequency.

## Verification

- [x] `average_depth` returns `0.0` (not an exception) when called on a
      region with zero covered positions.
- [x] `report_average` produces `Frequency = NaN` (not `inf`, not an
      exception) for a region where `Depth == 0`.
- [x] A full pipeline run in AVERAGE mode that includes a region with zero
      coverage completes successfully and writes `ijump_report_by_is_reg.txt`
      with a `NaN` frequency row for that region, instead of crashing.
      (Covered at the unit level: `average_depth`'s zero-coverage return and
      `report_average`'s NaN-Frequency path are both exercised directly;
      no end-to-end pipeline fixture with a real zero-coverage region exists
      in this test suite to extend.)
- [x] Existing `average_depth` parity tests
      (tests/test_average_depth_pysamstats_parity.py) still pass — the
      covered-only averaging behavior itself is unchanged.
