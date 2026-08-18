# 04 — Stop `report_average` pretending its caller controls the minimum match length

**What to build:** The average-mode report function's signature tells the truth
about what actually affects its output. Today it accepts a minimum-match-length
argument, and then throws it away — so a reader (or an agent) tuning that value
upstream would see no effect and have no idea why.

## Why

`report_average` in `region_summary.py` declares a `min_match` parameter, then
overwrites it on the first line of the body with `min_match = min(match_lengths)`
before any read. `ISClipped` sets `self.min_match = 150` and dutifully passes it
at the call site, where it is silently discarded.

The behaviour is correct and matches `master` — this is not a bug report. The
problem is that the seam lies: `min_match` looks like a knob, and isn't. The
review flagged it on both axes (as dead-parameter/Speculative Generality, and as
a misleading seam left by the ticket-13 extraction).

The same file has a smaller instance in `summarize_by_region`: a local `columns`
list is built with `columns = [...]` then `columns.extend(...)` and never read.
Ruff's F841 misses it precisely because the `.extend` call counts as a use, so
tooling will not catch this one.

**Note the adjacent trap:** the neighbouring `blast_min` parameter *is* genuinely
used, in the frequency calculation further down the same function. Do not remove
it by pattern-matching on "unused-looking numeric parameter".

## Scope

- Remove the `min_match` parameter from `report_average` and drop it from the
  call site.
- **Decide what happens to `ISClipped.min_match` and record the decision.** Two
  honest options: (a) it has no other consumer, so delete the attribute too; or
  (b) it is meant to be a real knob, in which case `report_average` should honour
  the passed value instead of recomputing it — but that **changes output**, so it
  becomes its own ticket, not a quiet edit inside this one. Grep for other
  consumers before choosing. Default to (a) unless a consumer turns up.
- Delete the dead `columns` local in `summarize_by_region`.

## Out of scope

- Touching `blast_min` — it is live.
- Any change to the computed frequency values. This ticket is strictly
  behaviour-preserving; if option (b) above looks right, file it separately.

## Verification

- `tests/test_region_summary.py` passes unchanged. It covers this function, so it
  is the safety net here — if a test needed editing, the change was not
  behaviour-preserving and something is wrong.
- Report output for a real sample is byte-identical before and after. Given this
  is meant to be a pure signature cleanup, that is the whole bar.
- `ruff`/`mypy` still pass clean on `src/ijump/`.
- The decision about `ISClipped.min_match` is written in Comments, whichever way
  it went.

**Blocked by:** None — can start immediately.

**Status:** ready-for-agent

- [ ] `min_match` parameter and its argument at the call site removed.
- [ ] `ISClipped.min_match` grepped for other consumers; decision made and recorded.
- [ ] Dead `columns` local removed from `summarize_by_region`.
- [ ] `blast_min` left alone.
- [ ] `tests/test_region_summary.py` passes with no edits.
- [ ] Report output byte-identical on a real sample.

## Comments
