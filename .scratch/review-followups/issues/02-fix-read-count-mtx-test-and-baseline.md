# 02 — Fix the test this branch broke, and retract the "pre-existing failure" claim

**What to build:** `pytest` passes from a clean clone with no known failures,
and no ticket file in `.scratch/` still tells a future reader that the failure
was pre-existing and unrelated. Right now one test fails, and the tracker
records that failure as somebody else's problem.

## Why

Two separate problems, one root cause.

**The test.** `tests/test_no_results_paths.py` has a test asserting that the
read-count-matrix helper rejects an invalid orientation argument. It calls that
helper as a static method **on `ISClipped`**. The isclipped-refactor ticket 09
extraction moved the helper out to module level in `frequency_estimation.py` and
left no delegating alias behind, so the call raises `AttributeError` where the
test asserts `pytest.raises(ValueError)`. Ticket 09's own "Done when" required
*"`pytest` passes from a clean clone"*; it does not.

**The false baseline — the more damaging half.** Several later tickets record
this failure in their Comments as a *"pre-existing, unrelated failure … references
a method that no longer exists"*. It is not pre-existing: `master` has the helper
as an `ISClipped` static method, and the test passes there. This branch broke it.

The claim then propagated. `.scratch/packaging/issues/08-upgrade-numpy-2.md`
bakes it into its **Verification** section as the expected baseline ("45 passed,
1 pre-existing unrelated failure in …"), so an agent following that ticket was
instructed to expect the failure and move on. Every subsequent ticket inherited
a test baseline that was wrong about its own provenance.

## Scope

- Point the test at the helper's real home in `frequency_estimation`. Confirm
  the assertion still tests what it was written to test — the invalid-orientation
  guard — and not merely that *some* exception is raised.
- Decide and record whether a private module-level function should be under test
  directly at all, or whether the test belongs against the public
  frequency-estimation entry point. Either answer is fine; state which and why.
- **Correct the record.** Find every ticket file referencing this failure —
  `grep -rln "pre-existing" .scratch/*/issues/` is the starting point (13 files
  match; the review confirmed 8 make the pre-existing claim, so read each rather
  than editing all 13). For each, **append a correction to `## Comments`**
  stating that the failure was introduced by isclipped-refactor ticket 09, not
  pre-existing, and that it is fixed by this ticket.
- **Do not rewrite the original ticket bodies.** Commits `1096b15` and `c506b75`
  ("Restore ticket NN's original text after merge, keep agent's findings")
  establish that the maintainer restores original ticket wording and keeps agent
  findings appended separately below. Honour that: corrections go in Comments.
  The one exception is packaging/08's **Verification** section, which states the
  wrong expected baseline as an instruction to future implementers — that line
  needs correcting in place, with a Comment noting the edit and why.
- Add a Comment to `.scratch/isclipped-refactor/issues/09-extract-frequency-estimation.md`
  recording that its "pytest passes from a clean clone" criterion was not
  actually met at the time it was closed.

## Out of scope

- Reverting or redesigning ticket 09's extraction. Moving the helper to module
  level was correct; only the un-updated caller was wrong.
- Re-running the full suite in a working conda environment, or unticking
  packaging/08's overstated `pytest` checkbox — that is ticket 08's job.

## Verification

- The previously failing test passes.
- `grep -rn "pre-existing" .scratch/*/issues/` surfaces no remaining claim that
  *this particular* failure was pre-existing. (Other, genuinely pre-existing
  issues may legitimately still use the phrase — check what each hit refers to.)
- No original ticket body text was altered except packaging/08's incorrect
  Verification baseline, and that edit is noted in its Comments.
- Full-suite result recorded honestly in Comments, including how many files
  could not be collected for want of `pysamstats` in the sandbox — that
  limitation is real and documented, and this ticket should not pretend
  otherwise.

**Blocked by:** None — can start immediately.

**Status:** ready-for-agent

- [ ] The test calls the helper at its real location and passes.
- [ ] A decision is recorded on whether to test the private helper or the public entry point.
- [ ] Every ticket claiming this failure was pre-existing carries a Comment retracting that.
- [ ] packaging/08's incorrect expected-baseline line in Verification is corrected in place and the edit noted.
- [ ] isclipped-refactor ticket 09 carries a Comment about its unmet "pytest passes" criterion.
- [ ] No original ticket body text rewritten beyond that one baseline line.

## Comments
