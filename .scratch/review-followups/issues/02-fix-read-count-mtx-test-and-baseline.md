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

**Status:** done

- [x] The test calls the helper at its real location and passes.
- [x] A decision is recorded on whether to test the private helper or the public entry point.
- [x] Every ticket claiming this failure was pre-existing carries a Comment retracting that.
- [x] packaging/08's incorrect expected-baseline line in Verification is corrected in place and the edit noted.
- [x] isclipped-refactor ticket 09 carries a Comment about its unmet "pytest passes" criterion.
- [x] No original ticket body text rewritten beyond that one baseline line.

## Comments

Implemented 2026-08-17.

**The test fix.** `tests/test_no_results_paths.py::test_read_count_mtx_rejects_invalid_orientation`
now imports and calls `frequency_estimation._read_count_mtx` directly instead
of the no-longer-existing `ISClipped._read_count_mtx`, and asserts on the
guard's actual message (`pytest.raises(ValueError, match='"left" or "right"')`)
so it keeps testing the invalid-orientation guard specifically, not just "some
exception."

**Decision: test the private helper directly, not the public entry point.**
`_read_count_mtx` has no public wrapper. `estimate_frequencies` (the module's
public entry point) always calls it with a hardcoded `"left"` or `"right"` and
never exposes `orientation` to its own callers — the invalid-orientation guard
is simply unreachable from the public API. Testing the private module-level
function directly is therefore the only way to exercise that guard at all, so
that's what the fixed test does (this mirrors what the pre-refactor
`ISClipped._read_count_mtx` static-method test already did — same shape of
test, just pointed at the helper's new home).

**Correcting the record.** `grep -rln "pre-existing" .scratch/*/issues/`
returned 13 files. Read each rather than editing blindly; 8 actually make the
specific claim that `test_read_count_mtx_rejects_invalid_orientation` was
pre-existing/unrelated:
- `.scratch/isclipped-refactor/issues/12-close-circos-av-depth-seam.md`
- `.scratch/isclipped-refactor/issues/16-evaluate-dropping-pysamstats.md`
- `.scratch/packaging/issues/01-src-layout-package.md`
- `.scratch/packaging/issues/02-cli-subcommand-dispatch.md`
- `.scratch/packaging/issues/03-uv-migration.md`
- `.scratch/packaging/issues/04-local-conda-packaging.md`
- `.scratch/packaging/issues/06-ruff-mypy-config.md`
- `.scratch/packaging/issues/08-upgrade-numpy-2.md`

Each got a Comment appended (body text otherwise untouched) retracting the
claim and pointing at this ticket as the fix. The other 5 matches were false
positives on the same grep string, checked and left alone because they're
about something else entirely:
- `.scratch/isclipped-refactor/issues/06-find-pair-characterization-test.md`
  — "pre-existing environment gap" = missing `pysamstats`/`Bio`/`sklearn` in
  the shell, unrelated to this test.
- `.scratch/isclipped-refactor/issues/09-extract-frequency-estimation.md`
  — "No pre-existing captured-real-data fixture" = about test-fixture
  provenance, not a failure claim. (It still got a *separate* Comment per
  this ticket's own scope, about its unmet "pytest passes from a clean
  clone" criterion — see below.)
- `.scratch/isclipped-refactor/issues/14-fix-estimation-mode-default.md`
  — "Confirmed pre-existing" refers to the `--estimation_mode` default bug,
  a different failure entirely.
- `.scratch/isclipped-refactor/issues/15-commit-untracked-agent-scaffolding.md`
  — same `pysamstats` environment-gap note as ticket 06, not about this test.
- `.scratch/packaging/issues/07-lint-pre-commit-ci.md` — "pre-existing
  findings" refers to ruff lint findings in `tests/`, unrelated.

`.scratch/packaging/issues/08-upgrade-numpy-2.md` was the one case where the
false claim sat in a **Verification** section as an instruction to future
implementers ("pytest passes with the same profile as today (45 passed, 1
pre-existing unrelated failure...)"), so that line was corrected in place (to
"46 passed, 0 failures") per this ticket's explicit exception, with a Comment
in that file noting the edit and why.

`.scratch/isclipped-refactor/issues/09-extract-frequency-estimation.md` got a
Comment recording that its "Done when" criterion "`pytest` passes from a
clean clone" was not actually met at close time — this is the ticket whose
extraction broke the test.

**Sandbox limitations, recorded honestly.** This sandbox's system Python
(3.9.13, Homebrew) has a broken/ABI-mismatched pandas/pyarrow install, so a
fresh venv was needed. In that venv, `pysam`/`pysamstats` could not be built:
`pysamstats`'s isolated build pulls in `pysam<0.16` from source, which fails
on `ModuleNotFoundError: No module named 'pkg_resources'` (newer setuptools
dropped it) even after pre-installing an older setuptools in the venv itself,
because `pysamstats`'s own build-isolation environment ignores that; forcing
`--no-build-isolation` instead hits `no cython installed, but can not find
pysam/libchtslib.c` since this pysam release ships no compiled sources named
that. This reproduces the same "conda-only, unbuildable via plain pip/venv"
wall packaging tickets 03/04 already documented. To get `tests/` collecting
at all, a **local stub** `pysamstats` module (raising `NotImplementedError`
from `load_coverage`) was placed in the venv's site-packages — not part of
this ticket's diff, not committed, exists only in `/tmp/ijump-venv-t02` for
this run.

With that stub, `pytest -q` from this venv: 36 passed, 10 failed. All 10
failures are attributable to the sandbox, not to this ticket's change or to
real bugs:
- `tests/test_average_depth_pysamstats_vs_count_coverage.py` (3 tests) and
  `tests/test_isclipped.py::test_average_depth_returns_mean_coverage_for_region`
  — hit the stub's `NotImplementedError` (real `pysamstats.load_coverage`
  unavailable).
- `tests/test_empty_run_outputs.py` (4 tests),
  `tests/test_estimation_mode_default.py::test_omitted_estimation_mode_writes_full_average_mode_output_set`,
  and `tests/test_estimation_mode_validation.py::test_invalid_estimation_mode_rejected_at_parse_time`
  — `subprocess`-invoke the `ijump` console script, which isn't installed
  because `pip install -e .` also depends on the same unbuildable
  `pysam`/`pysamstats` chain.

The specific test this ticket targets,
`tests/test_no_results_paths.py::test_read_count_mtx_rejects_invalid_orientation`,
passes, along with the rest of `tests/test_no_results_paths.py` (4/4) and the
full suite modulo the 10 sandbox-attributable failures above. Per this
ticket's own out-of-scope list, the full suite was not re-run in a working
conda environment and packaging/08's pytest checkbox was left untouched —
that remains ticket 08's job.

**Review note.** Two parallel sub-agent code-review tasks (Standards axis,
Spec axis) were launched against this diff per the `/implement` workflow.
Their `SendMessage` replies back to this session's agent name failed to
route; the coordinator relayed the Standards result directly instead. That
result: no hard violations — every `.scratch/*/issues/*.md` correction is a
pure Comments append except packaging/08's authorized in-place Verification
edit (which correctly carries an "Edit note:" disclosing it); 8 of the 13
`grep -rln "pre-existing"` hits are the genuine claim and were corrected, the
other 5 were correctly left alone; the test fix correctly repoints to
`frequency_estimation._read_count_mtx` with a tightened assertion. Sole
observation (not a violation): the same correction paragraph is duplicated
verbatim across the 7 append-only files — defensible given the repo's
append-only Comments convention.

The Spec-axis sub-agent also returned (relayed the same way): no missing
requirements, no scope creep, the checklist was correctly checked off, and
the test fix correctly targets the invalid-orientation guard. One note from
that review: re-running `grep -rln "pre-existing" .scratch/*/issues/` now
returns 14 files, not the ticket's stated 13 — the extra hit is this ticket
file itself (`review-followups/issues/02-...md`), whose own prose uses the
word "pre-existing" several times while describing the bug. Not an
implementer error; the ticket's "13 files" count was taken before this
ticket's own file existed in the count.

Manual self-review performed in parallel: re-read the full diff (`git diff
HEAD`) hunk by hunk, confirmed `tests/test_no_results_paths.py`'s fix
imports and calls `frequency_estimation._read_count_mtx` (verified against
the real current location via `grep -n "_read_count_mtx"
src/ijump/frequency_estimation.py`) and that the assertion still targets the
invalid-orientation guard specifically (`pytest.raises(ValueError,
match='"left" or "right"')`, not a bare `pytest.raises(Exception)`), re-ran
`tests/test_no_results_paths.py` (4/4 pass), and ran `ruff check` on the
changed test file — it flags one `I001` unsorted-import finding, but that
same finding already exists on `HEAD` before this change (verified via `git
show HEAD:... | ruff check --stdin-filename ... -`), so it's pre-existing
and, per ticket 07's Comments, `tests/` is out of ruff's currently-enforced
scope anyway — not something this ticket introduced or needs to fix.
Confirmed via `git diff --stat HEAD` that every `.scratch/*/issues/*.md`
change other than packaging/08's one authorized Verification line is a pure
append (no lines removed from original ticket bodies).

**Status corrected 2026-08-25.** The work above landed on 2026-08-17 with every box ticked,
but the `Status:` line was left at `ready-for-agent`, so the ticket kept showing up on the
agent frontier. Re-verified before relabelling: the test imports and calls
`frequency_estimation._read_count_mtx`, asserts on the guard's own message, and
`tests/test_no_results_paths.py` passes 4/4; `packaging/08`'s Verification baseline reads
"46 passed, 0 failures". Nothing to re-do.
