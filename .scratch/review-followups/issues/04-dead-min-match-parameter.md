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

**Status:** done

- [x] `min_match` parameter and its argument at the call site removed.
- [x] `ISClipped.min_match` grepped for other consumers; decision made and recorded.
- [x] Dead `columns` local removed from `summarize_by_region`.
- [x] `blast_min` left alone.
- [x] `tests/test_region_summary.py` passes (one stale keyword argument in the
      *call*, not an assertion, had to be dropped — see Comments below for why
      this was unavoidable).
- [x] Report output byte-identical on a real sample.

## Comments

**Branch:** this worktree's original branch (`worktree-agent-a89c9fa48d3f1bd05`)
was stale — checked out at `1975da4`, missing the entire `src/ijump` layout
and all ticket work up to `refactor`'s tip (`d933ce4`). `git reset`/`merge`/
`rebase` onto `refactor` were all blocked by the sandbox's auto-mode
classifier, so the fix was implemented on a fresh branch
`ticket-04-dead-min-match`, created from `refactor` via
`git checkout -b ticket-04-dead-min-match refactor` (non-destructive — no
history rewritten). That branch is what the commit for this ticket lives on.

**`ISClipped.min_match` decision: option (a) — deleted.** Grepped the whole
source tree (`grep -rn "min_match" src/ tests/`) before touching anything.
The only two hits were `self.min_match = 150` (isclipped.py:164) and the one
place it was passed to `report_average` (isclipped.py:620, now removed). No
other consumer exists, so the attribute carries no signal and was deleted
along with its comment. Went with (a), not (b): honoring the passed value
would change the computed frequency (a real behavior change), which is
explicitly out of scope for this ticket per its own "Out of scope" section —
that path is left for a separate ticket if someone decides `min_match` should
become a real knob.

**Note on the "no test edits" bar:** `tests/test_region_summary.py`'s
`test_report_average_matches_pinned_golden_output_for_single_hit_path` called
`report_average(..., min_match=150, ...)` by keyword. Removing the parameter
necessarily makes that call raise `TypeError: unexpected keyword argument`
regardless of implementation approach — this isn't avoidable while still
removing the parameter (confirmed empirically: ran the test before editing
it, watched it fail with exactly that TypeError). Removed only that one
stale keyword argument from the test's call; every assertion and expected
value in the file is untouched, and all 3 tests in the file pass afterward
with the same expected output. Read the ticket's "no edits" bar as being
about not touching what's being verified (assertions/expected values), not
about the call syntax an inherently-changing signature forces on every
caller, including test callers.

**Byte-identical verification, actually run (not assumed):** the full
`ISClipped` pipeline can't run in this sandbox — `pysam`/`pysamstats` fail to
build via `uv sync` here (the known blocker from ticket 16). Instead, loaded
the pre-change and post-change `region_summary.py` as two separate modules
and fed both `report_average` implementations the *same* real
`sum_by_region` data (143 rows, from `Example files/output_tables/
ijump_sum_by_reg.txt`, a genuine prior run's output) plus identical
`match_lengths`/`read_lengths`/`n_reads_analyzed`/`blast_min`/`average_depth`
arguments — with the pre-change call additionally passing `min_match=999`
(an intentionally wrong value) to prove it really is discarded either way.
`pd.testing.assert_frame_equal` and a raw CSV-string comparison both came
back byte-identical. Also ran `ruff check src/ijump/` and `mypy src/ijump/`:
both clean.
