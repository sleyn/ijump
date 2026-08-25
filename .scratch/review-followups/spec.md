# Review follow-ups

Follow-up work from a two-axis code review (Standards + Spec) of the `refactor`
branch against `master` — 40 commits, 91 files, +8507/-2428 — run 2026-08-17.
That range covers both earlier ticket batches: `.scratch/isclipped-refactor/`
(tickets 01-16) and `.scratch/packaging/` (tickets 01-08).

These are follow-ups, not new features. Nothing here is a request to re-litigate
a decision either batch already made.

## Priority

Two findings are genuine bugs that ship broken behaviour today:

- **01** — `DataFrame.append` was a latent upgrade blocker while pandas was
  pinned to `1.3.5`. The packaging batch widened the pin to `pandas<3`, which
  turns it into a live `AttributeError` in a wired-up CLI subcommand.
- **02** — a test broken by this branch's own ticket-09 extraction was recorded
  as a "pre-existing, unrelated failure" in the Comments of several later
  tickets. `master` passes it. A false baseline is now embedded in the tracker
  and was trusted by subsequent tickets.

Everything else is a tracker-hygiene gap or a code smell worth capturing before
it calcifies.

## Ticket numbering

**Reference these as `(review-followups NN)` in commit messages — never bare
`(ticket NN)`.** Numbers `01-08` already exist in both
`.scratch/isclipped-refactor/issues/` and `.scratch/packaging/issues/`, and the
review found that collision actively caused commits to be mis-mapped to the
wrong spec. Don't add a third ambiguous namespace.

## Tickets

| NN | Summary | Blocked by | Status |
| -- | ------- | ---------- | ------ |
| 01 | `DataFrame.append` breaks `isfinder-parse` under pandas 2 | — | wontfix (module deleted by `isfinder-annotation 02`) |
| 02 | Fix `test_no_results_paths` and correct the false baseline | — | ready-for-agent |
| 03 | Clear `tests/` lint findings, widen the pre-commit/CI gates | — | ready-for-agent |
| 04 | Drop the dead `min_match` parameter and dead local | — | ready-for-agent |
| 05 | Collapse two duplicated definitions into one home each | — | ready-for-agent |
| 06 | Naming-consistency sweep (typos, unfinished rename, glossary) | — | ready-for-agent |
| 07 | Move the Circos call assembly onto `ISClipped` | 06 | ready-for-agent |
| 08 | Conda/Docker verification debt and tracker status hygiene | 01, 02 | ready-for-human |
| 09 | Give the boundary clump a type; shrink a 13-parameter signature | 05, 07 | needs-triage |

01-06 are independent and can run in parallel.

**Two are deliberately not AFK-ready:**

- **08** needs a human with a working conda and Docker installation plus a real
  sample. Every agent session so far has hit the same wall (no conda in the
  sandbox; `pysamstats` cannot be pip-built — see `.scratch/packaging/issues/03-uv-migration.md`).
- **09** is a design proposal, not a defect. It needs a maintainer's scope call,
  and `.scratch/isclipped-refactor/notes/no-test-safety-net.md` documents that
  the repo lacks the coverage to make a change of that shape safely.

## Conventions

Line numbers in these tickets were verified against `HEAD` at filing time but
will drift. Re-grep for the named symbol rather than trusting the line —
`.scratch/packaging/issues/08-upgrade-numpy-2.md` sets that precedent.

Where a ticket finds a genuine bug outside its own scope, flag it in `## Comments`
rather than fixing it inline. That convention comes from
`.scratch/packaging/issues/06-ruff-mypy-config.md`, and the review found one
place where it was broken (see ticket 08).
