# 07 — pre-commit hook + lint-only CI workflow

**What to build:** Actual enforcement for ticket 06's ruff/mypy config: a
`.pre-commit-config.yaml` running both locally at commit time, plus a new
`.github/workflows/lint.yml` running both on every push/PR — the repo's
first CI workflow of any kind.

## Why

Grilled directly with the user. A local-only pre-commit hook can always be
bypassed (`--no-verify`, or a contributor who never ran `pre-commit
install`), so both were wanted together rather than pre-commit alone.

`pytest` was explicitly kept out of this workflow: this repo's tests need
BLAST+ and `pysamstats` (a compiled, conda-only dependency per
`README.md` — confirmed by reproducing `ModuleNotFoundError: No module
named 'pysamstats'` when running plain `pytest` in this environment
without conda), and there's no known-good CI install recipe for those yet
— that's what tickets 03 (uv)/04 (conda)/05 (Docker) are for. Bundling a
from-scratch BLAST+/pysamstats install into this lint ticket would
duplicate work those tickets are already scoped to do properly. Adding
`pytest` to CI is deferred to a future ticket once one of those lands and
has a reusable dependency-install recipe.

## Scope

- Depends on ticket 06 (the ruff/mypy config this ticket wires up must
  exist first).
- `.pre-commit-config.yaml` at the repo root: hooks for `ruff check --fix`,
  `ruff format`, and `mypy` (using each tool's official pre-commit-hooks
  repo mirror, e.g. `astral-sh/ruff-pre-commit`,
  `pre-commit/mirrors-mypy` — pin versions matching whatever ticket 06
  installed).
- Document `uv run pre-commit install` (or `pip install pre-commit &&
  pre-commit install` if uv isn't fully wired by the time this lands — use
  whichever install path ticket 03 has established) in the README's
  development section, so contributors know to run it once per clone.
- `.github/workflows/lint.yml`: triggers on push/PR, sets up Python
  (matching `pyproject.toml`'s `requires-python`), installs the package
  with dev/lint extras (however ticket 01/03 expose them — e.g. `uv sync
  --extra lint` or `pip install .[lint]`, whichever dependency-group
  mechanism those tickets chose), then runs `ruff check`, `ruff format
  --check`, `mypy src/ijump` as separate steps (so a failure in one is
  clearly attributed rather than buried in a single combined step).
- Explicitly no `pytest` step in this workflow — see Why. Leave a comment
  in the workflow file noting this is intentional and pointing at the
  ticket that should add it later, so it doesn't read as an oversight.

## Out of scope

- Adding `pytest` to CI (deferred, see Why — needs a BLAST+/pysamstats
  install recipe from tickets 03/04/05 first).
- Any other CI workflow (build/release automation, Docker image build
  checks, etc.) — this ticket is lint-only.
- Branch protection rules requiring this workflow to pass before merge —
  that's a repo-settings change outside this ticket's file-based scope;
  flag it as a follow-up if relevant once this lands.

## Verification

- `pre-commit run --all-files` passes clean locally after `pre-commit
  install`.
- A deliberately-introduced lint violation (e.g. an unused import) in a
  local test commit is caught and blocked by the pre-commit hook before
  the commit completes — confirm interactively, then revert the test
  violation.
- Pushing a branch/opening a PR triggers `.github/workflows/lint.yml` and
  it passes on the current (post-ticket-06) codebase.
- The same deliberately-introduced violation, pushed with `--no-verify` to
  bypass the local hook, is still caught by the CI workflow.

**Blocked by:** 06 (ruff + mypy config)

**Status:** done, mostly — everything file-based is done and verified locally; the
actual live push/PR trigger of `.github/workflows/lint.yml` on GitHub's runners
could not be confirmed from this sandbox (see Comments).

- [x] `.pre-commit-config.yaml` added, running ruff check/format + mypy.
- [x] README documents the one-time `pre-commit install` step.
- [x] `.github/workflows/lint.yml` added, running ruff check/format + mypy as separate steps.
- [x] No `pytest` step in this workflow; the omission is documented inline as intentional.
- [x] `pre-commit run --all-files` passes clean.
- [ ] CI workflow passes on the current codebase and is confirmed to catch a deliberately-introduced violation. **Not confirmed on real GitHub Actions runners** — this sandbox has no way to trigger an actual push/PR workflow run. Validated by inspection instead (see Comments); a maintainer should confirm on a real push/PR before relying on this as a merge gate.

## Comments

Implemented 2026-08-17, working from a worktree branched off `refactor` at
`4db4ee3` (ticket 06's tip). Environment: the same conda env ticket 06 verified
against, `/Users/sleyn/miniconda3_envs_backup/ijump-verify` (Python 3.11, ruff
0.16.3, mypy 2.3.1), plus `pre-commit` 4.6.2 installed into that env for this
ticket (not a runtime dependency, dev-only like ruff/mypy).

**Files added/changed:**

- `.pre-commit-config.yaml` — `astral-sh/ruff-pre-commit` `rev: v0.16.3` (`ruff`
  with `--fix`, `ruff-format`) and `pre-commit/mirrors-mypy` `rev: v2.3.1`
  (`mypy`, `args: [src/ijump]`), pinned to the exact versions ticket 06
  installed. All three hooks are scoped to `files: ^src/ijump/` — ticket 06's
  clean baseline only covers `src/ijump/`; `tests/` was explicitly out of that
  ticket's scope and is not currently ruff-clean (verified: `ruff check .` from
  repo root finds 17 pre-existing findings in `tests/`, none introduced by this
  ticket). Scoping the hooks avoids failing commits/CI on out-of-scope,
  pre-existing findings while still enforcing the real baseline.
- `pyproject.toml` — added `[project.optional-dependencies] lint = ["ruff==0.16.3", "mypy==2.3.1"]`.
  No such extra/dependency-group existed yet from tickets 01/03, so this ticket
  defines it, pinned to match `.pre-commit-config.yaml`'s versions.
- `.github/workflows/lint.yml` — new workflow, triggers on `push`/`pull_request`,
  Python 3.11 via `actions/setup-python`, then `ruff check`, `ruff format
  --check`, `mypy src/ijump` as three separate steps (all scoped to
  `src/ijump`, same reasoning as above). No `pytest` step, with an inline
  comment explaining why and pointing at tickets 03/04/05.
- `README.md` — new "Lint / pre-commit" subsection under Development setup,
  documenting `uv run pre-commit install` (with a `pip install pre-commit &&
  pre-commit install` fallback, since ticket 03 already documents that `uv
  sync` doesn't work end-to-end on this machine) and `pre-commit run
  --all-files`.

**Judgment calls / deviations:**

- **Python version in CI (3.11, not `requires-python`'s `>=3.7` floor).**
  `requires-python = ">=3.7"` is stale (ticket 06's own comments note it) —
  mypy 2.3.1 refuses `python_version` below 3.10, and `environment.yml`
  actually targets Python 3.11. Used 3.11 in `actions/setup-python` to match
  the environment the ticket 06 baseline was verified against, with an inline
  comment explaining why, same spirit as ticket 06 leaving `python_version`
  unset in `[tool.mypy]`.
- **CI installs `ruff==0.16.3 mypy==2.3.1` directly rather than `pip install
  '.[lint]'`.** Tested this explicitly: `pip install '.[lint]'` would also try
  to build the project's core runtime deps (`pysam`, `pysamstats`), which
  ticket 03 already found cannot be resolved/built via `pip` or `uv` on a
  current Python (unmaintained `pysamstats==1.1.2` hard-pins an unbuildable
  `pysam==0.15.4`). Confirmed neither `ruff` nor `mypy` need `ijump`'s runtime
  deps installed to lint/type-check the source tree: ran both tools from a
  brand-new venv with *only* `ruff==0.16.3`/`mypy==2.3.1` installed (no
  pandas/pysam/numpy/etc. present at all) against `src/ijump`, and both passed
  clean identically to the full conda env. Documented this reasoning inline in
  the workflow file so it doesn't look like an oversight of the `[lint]` extra
  added to `pyproject.toml`.
- **Scoped ruff/mypy to `src/ijump/` in both the pre-commit hooks and CI**,
  rather than the repo root the ticket's prose examples show (`ruff check`,
  `ruff format --check`). Necessary because `tests/` was never brought to a
  clean baseline by ticket 06 (confirmed 17 findings there, all pre-existing —
  unsorted imports, `E501` long lines, 2 unused imports); running unscoped
  would make both the pre-commit hook and CI red on day one for reasons
  outside this ticket's scope. `mypy`'s scope was already `src/ijump` per the
  ticket text; extended the same scoping to the two `ruff` steps/hooks for
  consistency. Flagging `tests/`'s lint cleanliness as a reasonable follow-up
  ticket if the project wants both directories enforced.

**Verification performed:**

- `ruff check src/ijump`, `ruff format --check src/ijump`, `mypy src/ijump` —
  all pass clean standalone (both in the `ijump-verify` conda env and in a
  from-scratch venv with only `ruff`/`mypy` installed, no runtime deps).
- `pre-commit install` — installs the git hook successfully (hooks live at the
  primary repo's `.git/hooks/`, shared across worktrees, since this worktree's
  `.git` is a worktree pointer file).
- `pre-commit run --all-files` — `ruff check`, `ruff format`, `mypy` all report
  `Passed`.
- Deliberate-violation check: appended a stray `import os` to
  `src/ijump/cli.py`, staged it, and ran `git commit`. The `ruff check` hook
  failed the commit (`Found 1 error (1 fixed, 0 remaining)` — auto-fixed and
  blocked so the fix can be reviewed rather than silently committed), which
  blocked the commit from completing. Reverted the test change afterward
  (`git checkout HEAD -- src/ijump/cli.py`); it is not part of this ticket's
  commit.
- `.pre-commit-config.yaml` / `.github/workflows/lint.yml` YAML parse cleanly
  (`yaml.safe_load`); `pre-commit validate-config` passes.
- **Not done / cannot be done from this sandbox:** actually pushing a branch or
  opening a PR to confirm `.github/workflows/lint.yml` runs and passes on
  real GitHub Actions infrastructure, or confirming it catches a `--no-verify`
  bypassed violation there. This needs a live push/PR against the GitHub repo,
  which this environment cannot perform. Everything short of that (config
  validity, the exact commands GitHub Actions would run, clean-pass behavior
  against the current codebase) was verified locally as closely as possible.
  A maintainer should confirm the actual Actions run once this branch is
  pushed.
- Explicitly did not add a `pytest` step, per scope — confirmed the workflow
  file only runs `ruff check` / `ruff format --check` / `mypy src/ijump`, no
  test invocation.
- Did not add branch protection rules or any other CI workflow — out of
  scope, untouched.
