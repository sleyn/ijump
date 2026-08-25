# 08 — Close the conda/Docker verification debt and make the tracker tell the truth

**What to build:** The packaging batch's claims are all actually true. Someone
picking up this repo can create the documented conda environment, build the conda
package, build and run the Docker image, and run a real sample through — and the
ticket files accurately say which of those has been verified and which hasn't.
Right now several boxes are ticked on work that no session was able to run, and
two implemented tickets still advertise themselves as available.

## Why

**This needs a human.** Every agent session on this batch hit the same wall and
documented it honestly: there is no working conda in the sandbox, and
`pysamstats` cannot be pip- or uv-installed because it hard-pins an unbuildable
`pysam` (see `.scratch/packaging/issues/03-uv-migration.md`, and packaging/08's
Comments). Agents substituted ad hoc venvs and said so. What is left genuinely
requires a machine with conda and Docker plus a real sample — it cannot be
delegated AFK, which is why this is `ready-for-human` and not `ready-for-agent`.

Three distinct things are wrong:

**Unrun verification.** packaging/08 and packaging/05 each leave real
checkboxes unticked — the numpy≥2 conda run and `conda build .` for 08, the
Docker equivalents for 05. Those are honest, and they are the actual work here.

**One overstated tick.** packaging/08 ticks *"pytest passes with no np.int0
deprecation warning"*. Its own Comments disclose that 6 of 14 test files could
not even be collected and only ~22 tests ran. The claim is broader than what was
verified.

**Status labels that mislead.** packaging/05 and packaging/08 are still
`Status: ready-for-agent` although both are implemented and merged — so the
frontier scan in `docs/agents/issue-tracker.md` would hand an agent work that is
already done. Separately, `done`, `done, partially`, and `done, mostly` are in
use across packaging 01-07 but appear nowhere in
`docs/agents/triage-labels.md`, whose table defines exactly five labels.

## Scope

- Create a real conda environment from `environment.yml`, confirm the solver
  actually resolves numpy ≥ 2 (don't infer it from the pin), and run the **full**
  suite there — the first time all 14 test files will have been collectable.
  Record the real numbers.
- Run `conda build .` against `meta.yaml`.
- Build the Docker image and run a real sample through it, per packaging/05's
  criteria.
- Run the manual real-sample `ijump run` end to end in the conda environment.
- Correct packaging/08's overstated `pytest` tick to match what is now actually
  true, and tick the boxes this work legitimately closes.
- Set `Status:` correctly on packaging/05 and packaging/08.
- **Resolve the label vocabulary**: either add the `done` variants to the table
  in `docs/agents/triage-labels.md` (it invites this — *"Edit the right-hand
  column to match whatever vocabulary you actually use"*), or migrate the
  packaging tickets to the five documented labels. Pick one; leave the tracker
  and its documentation agreeing.
- **Record one process deviation.** Commit `b6aca82` fixed a genuine pandas bug
  in `region_summary.py` under packaging/04, whose Scope listed only
  `environment.yml`, `meta.yaml`, and README. Ticket 06 states the convention:
  *"if something looks like a genuine bug … flag it in Comments rather than
  silently fixing it inline"* — and packaging/04's own Comments applies that
  convention to a different finding in the same breath, then breaks it for this
  one. **The fix is correct and covered by `tests/test_region_summary.py` —
  nothing to revert.** Just append a Comment to packaging/04 noting the
  deviation, so the record is accurate.

- **Run `pre-commit run --all-files`, and confirm CI is green.** Added 2026-08-25.
  `pre-commit` is not installed in any agent environment used on this repo, and no CI run is
  observable from one, so tickets 03 and 10 leave those two boxes unticked with the reason
  inline rather than ticked on trust. What *was* verified there: each hook's command
  (`ruff check`, `ruff format --check`, `mypy src/ijump tests`) passes on its own, and the
  shared `files: ^(src/ijump|tests)/` pattern matches test paths and not `simulation/`,
  `rule-tests/` or repo-root files. What is left is the hook runner itself and CI.

## Out of scope

- Reverting the `region_summary` fix.
- Fixing the broken test or the false "pre-existing failure" baseline — ticket 02
  owns that, and this ticket records the corrected result rather than deriving it.
- Bioconda submission. packaging/04 explicitly deferred that (weeks of upstream
  CI and maintainer review); this is local build verification only.

## Verification

- A conda environment created from `environment.yml` resolves numpy ≥ 2, and the
  full suite runs there with all test files collected. Actual pass/fail numbers
  recorded in Comments.
- `conda build .` succeeds.
- The Docker image builds and runs a real sample.
- A manual real-sample `ijump run` in the conda environment exits 0 and writes
  the expected output set.
- Every remaining ticked box in packaging/05 and /08 corresponds to something
  actually performed — re-read both checklists and confirm each one.
- `grep -rn "^\*\*Status:\*\*" .scratch/*/issues/` yields only labels present in
  `docs/agents/triage-labels.md`.

**Blocked by:** 02 (01 is closed `wontfix`).

- **01** — no longer a blocker: `isfinder-annotation 02` deleted the
  `isfinder-parse` subcommand entirely, so there is no crashing subcommand left to
  certify around and 01 is closed `wontfix`.
- **02** — the corrected pytest baseline is what this ticket records. Doing it
  first would mean writing down the wrong numbers again.

**Status:** ready-for-human

- [ ] `pre-commit run --all-files` passes with the hooks installed (carried from tickets 03 and 10).
- [ ] CI observed green (carried from ticket 03).
- [ ] Conda environment created; numpy ≥ 2 confirmed resolved, not assumed.
- [ ] Full test suite run in that environment with all 14 files collected; real numbers recorded.
- [ ] `conda build .` succeeds.
- [ ] Docker image builds and runs a real sample.
- [ ] Manual real-sample `ijump run` verified in the conda environment.
- [ ] packaging/08's overstated `pytest` tick corrected.
- [ ] `Status:` corrected on packaging/05 and packaging/08.
- [ ] Label vocabulary reconciled between the tracker and `triage-labels.md`.
- [ ] packaging/04 carries a Comment recording the out-of-scope bug fix.

## Comments
