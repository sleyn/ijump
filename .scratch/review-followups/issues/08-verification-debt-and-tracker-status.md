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

**Status:** done

- [x] `pre-commit run --all-files` passes with the hooks installed (carried from tickets 03 and 10).
- [x] CI observed green (carried from ticket 03) — `refactor` pushed 2026-08-27; Lint workflow run `33076602903` succeeded against `d9327bf`, the tip, matching `origin/refactor`.
- [x] Conda environment created; numpy ≥ 2 confirmed resolved, not assumed.
- [x] Full test suite run in that environment with all 14 files collected; real numbers recorded. (The repo now has 29 test files, not 14 — see Comments for the up-to-date count.)
- [x] `conda build .` succeeds.
- [x] Docker image builds and runs a real sample.
- [x] Manual real-sample `ijump run` verified in the conda environment.
- [x] packaging/08's overstated `pytest` tick corrected.
- [x] `Status:` corrected on packaging/05 and packaging/08.
- [x] Label vocabulary reconciled between the tracker and `triage-labels.md`.
- [x] packaging/04 carries a Comment recording the out-of-scope bug fix.

## Comments

**Closed out 2026-08-26, real conda + Docker available this session** (unlike
every prior agent session on this batch, which correctly reported no working
conda/Docker in their sandboxes — this one had both, so the "needs a human"
gate was actually reachable by an agent this time):

- `conda env create -f environment.yml` resolves; `numpy` lands at 2.x
  (`2.5.2` via `conda run`, `2.4.6` via the env's own interpreter directly —
  both ≥ 2, confirming the solver actually picked a numpy 2 release rather
  than an unpinned line happening to land on 1.x).
- Full suite (`pytest -m "not e2e"`) in that env: **194 passed, 4 skipped, 8
  deselected**, all 29 `tests/test_*.py` files collected (`198/206
  collected`). The repo has grown since this ticket was filed — it's 29 test
  files now, not 14 — but the point (every file collectable, none silently
  skipped by a broken import) holds and is stronger than before: the 6 files
  that couldn't be collected in prior ad hoc venvs needed `pysamstats`, which
  is no longer a dependency at all (isclipped-refactor ticket 16 round 2,
  `80f4235`).
- `pytest -m e2e` — the real-sample tier, `Test/Sample.bam` (840 MB) +
  `A_baumannii_assembly.fna` + GFF + IS table — **8 passed**, both estimation
  modes, byte-identical to committed goldens. This is what closes both the
  "manual real-sample `ijump run`" and (together with the plain-suite run
  above) the "full test suite" boxes: `golden_support.run_e2e_pipeline`
  actually shells out to the `ijump` CLI and asserts `returncode == 0`, so
  this is a real CLI invocation, not just library-level testing.
- `conda build .` succeeds against `meta.yaml`, producing
  `ijump-1.0.4-py_0.conda`; the recipe's own `ijump --help` test step passes.
- `docker build --platform linux/amd64 -t ijump:release-verify .` succeeds
  from the current tree; `docker run` against `tests/fixtures/` (mounted
  read-only, output mounted read-write) exits 0 and writes the expected
  output set.
- `pre-commit run --all-files` passes (`ruff check`, `ruff format`, `mypy`
  all green). Needed a throwaway Python 3.10+ venv for the
  `mirrors-mypy`/`ruff-pre-commit` hook environments — the machine's default
  `python3` is 3.9, and `mypy==2.3.1` (the pin in `.pre-commit-config.yaml`)
  requires ≥ 3.10. `python3.12` was available via Homebrew; installing
  `pre-commit` into a venv built on that and pointing it at the repo (rather
  than at the system Python) is what made the hook environments installable.
  This is a one-time local machine-setup detail, not a repo change.
- packaging/08's overstated `pytest` tick, packaging/05 and packaging/08's
  stale `Status: ready-for-agent`, the label-vocabulary table, and
  packaging/04's missing process-deviation comment are all corrected in
  those tickets directly (see their own files).

**Genuine finding surfaced during Docker re-verification, not fixed inline:**
`Dockerfile`'s run-dependency install (`biopython==1.79`, `"numpy<2"`,
`pysam==0.15.4`, plus the `build-essential`/`zlib1g-dev` apt packages) is
stale relative to current `src/ijump` — `pysamstats` (the reason for the
`pysam==0.15.4`/Python-3.8/apt-toolchain machinery) and the biopython BLAST
wrapper (the reason for pinning biopython) were both dropped from the
codebase after packaging/05 was written. The image still builds and runs
correctly, since carrying unused pinned packages isn't itself broken, but the
whole scaffolding is now removable complexity. Recorded in packaging/05's
Comments; left as a separate future ticket rather than fixed here, since this
ticket's own scope is verification, not re-architecture.

**Still open: CI observed green.** `git status` shows the `refactor` branch 4
commits ahead of `origin/refactor` at session start (`59d55eb` is the local
tip; `8b6063e`, the last pushed commit, is what CI last ran against — green).
Pushing those 4 commits and confirming CI on the actual tip is a `git push`,
which this session left for a human per this repo's convention of not taking
actions with effects beyond the local checkout without asking first. That's
the one box release-2.0.0/01 should still treat as open before tagging.
