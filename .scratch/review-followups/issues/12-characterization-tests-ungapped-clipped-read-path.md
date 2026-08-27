# 12 — Characterization tests over the ungapped clipped-read search path

**What to build:** `_crtable_ungapped` in `clipped_read_search.py` gets the same
kind of characterization coverage that junction pairing and region summary
already have from the isclipped-refactor batch — tests that pin its current
input/output behaviour (including the in-place mutation of its parameters)
before anything about its signature changes.

**Why:** `.scratch/isclipped-refactor/notes/no-test-safety-net.md` documents
that this path has no coverage, and review-followups/09 (now split into this
ticket plus 13 and 14) named that gap as the reason a 13-parameter,
in-place-mutating function couldn't be safely refactored yet. This mirrors how
isclipped-refactor tickets 01 and 06 preceded their own extractions.

**Blocked by:** None — can start immediately (09's own blockers, tickets 05 and
07, are both already done).

- [x] Characterization tests exist for `_crtable_ungapped` covering its normal
      path, its boundary/edge cases, and the in-place mutation of the
      parameters it currently mutates
- [x] Tests pass against the current, unmodified implementation
- [x] `ruff`/`mypy` clean on the new test file
- [x] Test file lives alongside the other `clipped_read_search.py` tests,
      following their existing naming/fixture conventions (added to
      `tests/test_clipped_read_search.py`, its own section)

**Unblocks:** review-followups/13 (the Boundary-type refactor depends on this
coverage existing first).

## Comments

**Verified against the original signature before 13/14 touched anything**
(same session): the tests below were written and run green against the
unmodified `_crtable_ungapped` first (`pytest tests/test_clipped_read_search.py`
-- 21 passed), confirming the pinned values before any refactor. 13 then
adapted the boundary construction in these same tests from positional
`[start, stop, edge, is_name, chrom]` lists to `Boundary(...)` calls (required
by 13's own signature change, not a behavior change -- full suite re-ran green
at 203 passed plus the e2e byte-identical check at that checkpoint), and 14
further rewrote the mutation-specific assertions into return-value assertions
plus one new non-mutation test, exactly as 14's own checklist anticipated
("extended if needed to assert non-mutation"). The tree's final state is the
post-14 version of these tests; the pre-refactor mutating contract they
originally pinned was verified interactively during implementation rather than
left frozen in the committed diff.
