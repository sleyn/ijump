# 03 — Bring `tests/` under the lint and format gates

**What to build:** A contributor who writes a badly formatted or import-messy
test gets told about it — by the pre-commit hook locally and by CI on push.
Today `tests/` is invisible to both, so lint debt accumulates there freely.

## Why

Packaging ticket 07 stated its acceptance criteria as CI running `ruff check`,
`ruff format --check`, and `mypy src/ijump`. What landed scopes **every**
pre-commit hook to `files: ^src/ijump/` and runs each CI command against
`src/ijump` only.

The narrowing is deliberate and documented in comments in both
`.pre-commit-config.yaml` and `.github/workflows/lint.yml`: ticket 06
established its clean baseline for `src/ijump/` alone, and gating on `tests/`
would have failed on pre-existing findings. That was a reasonable call at the
time, but the effect is permanent — nothing will ever pull `tests/` clean,
because nothing ever checks it.

The debt is small and worth clearing now rather than at 10x the size:
`ruff check tests/` reports **17 findings, 8 autofixable**, and they break down
as **9 × E501** (line-too-long), **6 × I001** (unsorted imports), **2 × F401**
(unused import). The 9 manual ones are pure line-wrapping. This is an afternoon,
not a project.

## Scope

- Clear the 17 findings in `tests/`. Autofix the import ones; wrap the long
  lines by hand. Per the convention in ticket 06, **list any non-autofix
  resolution in Comments** so a reviewer can check the judgement calls rather
  than re-deriving them from the diff.
- Widen the `files:` patterns on the `ruff` and `ruff-format` pre-commit hooks,
  and the corresponding CI commands, to cover `tests/` as well as `src/ijump/`.
  Replace the now-inaccurate explanatory comments in both files rather than
  leaving them contradicting the config.
- **Make an explicit decision about `mypy` and state it.** The mypy hook passes
  `args: [src/ijump]` with `pass_filenames: false`, so widening it is a separate
  choice from widening ruff. Test files are often the worst case for type
  checking (fixtures, monkeypatching, deliberately wrong inputs). Either bring
  `tests/` in and clear whatever it surfaces, or leave mypy scoped to
  `src/ijump/` **with the reason written in the config comment**. Do not leave
  it ambiguous.
- Re-check the finding counts before starting — they will have drifted.

## Out of scope

- Adding a `pytest` step to CI. `lint.yml` already documents why that is
  deferred: the suite needs BLAST+ and `pysamstats`, and no known-good CI
  install recipe exists yet — that waits on tickets 03/04/05 of the packaging
  batch.
- Other unlinted directories (`simulation/`, `rule-tests/`, repo-root scripts).
  If any are worth gating, note them in Comments as a follow-up; don't expand
  this ticket.
- Changing which ruff rules are enabled. This ticket applies the existing
  configuration to more files; it does not retune it.

## Verification

- `ruff check tests/` and `ruff format --check tests/` both pass clean.
- The chosen mypy scope passes, and the config comment states which scope was
  chosen and why.
- `pre-commit run --all-files` passes.
- A deliberately malformed test file is caught by the hook — confirm the widened
  `files:` pattern actually matches, rather than assuming the regex is right.
- CI passes on push.

**Blocked by:** None — can start immediately. (Ticket 02 also edits a file under
`tests/`; if both are in flight, expect a trivial conflict in that one file —
neither gates the other.)

**Status:** ready-for-agent

- [ ] All ruff findings in `tests/` resolved; non-autofix ones listed in Comments.
- [ ] `ruff` and `ruff-format` hooks and CI commands cover `tests/`.
- [ ] Stale "scoped to src/ijump because tests/ isn't clean" comments removed or rewritten in both config files.
- [ ] An explicit, documented decision made on mypy's scope.
- [ ] `pre-commit run --all-files` passes.
- [ ] Widened `files:` pattern confirmed to actually match a test file.
- [ ] CI green.

## Comments
