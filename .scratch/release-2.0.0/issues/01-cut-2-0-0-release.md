# 01 — Cut the 2.0.0 release

**What to build:** iJump 2.0.0 exists as a real, tagged release: `refactor`
(90 commits ahead of `master`, holding the entire IS-annotation subsystem, the
`ISClipped` module split, and the batch of bugfixes/audits since) is merged
into `master`, `pyproject.toml`'s version is bumped to `2.0.0`, `CHANGELOG.md`'s
`## Unreleased` section is retitled `## 2.0.0` with the Circos removal added,
and `master` is tagged `v2.0.0`.

**Why:** `pyproject.toml` already reads `1.0.4` but no tag or changelog entry
exists for it, and none of `refactor`'s work — including three already-known
breaking changes plus the Circos removal decided in this planning session —
has ever reached `master`. Semver calls for a major bump given the breaking
changes.

**This needs a human.** Merging 90 commits, running real conda/Docker
verification, and pushing a tag are exactly the kind of hard-to-reverse,
machine-dependent actions this repo's conventions reserve for a person, not an
AFK agent.

**Blocked by:**
- `review-followups/08` — the conda/Docker verification-debt ticket. Its full
  checklist (fresh conda env + full test suite, `conda build .`, Docker
  build+run, manual real-sample `ijump run`, tracker cleanup) is a hard gate:
  do not tag before it's clean.
- `circos-removal/01` — Circos must actually be gone before the release diff
  is final.
- `average-depth-zero-coverage/02` — its two remaining open questions
  (`blast_min`/`av_read_len` bias; interaction with `Depth=0`/`NaN`) must be
  resolved first, per the release decision to not defer them.
- The CHANGELOG draft (done directly in this session, not a separate ticket).

**Not blocked by:** `review-followups/12`, `13`, `14` (the Boundary-type
refactor and mutation removal split off from ticket 09) — that's a
design/readability improvement on working code, explicitly decided to land on
`master` after the tag rather than gate it.

**Status:** ready-for-human

- [x] `review-followups/08`'s full checklist passes — real conda/Docker verification
      done 2026-08-26/27: env, full suite (194 passed), `conda build .`, Docker
      build+run, manual real-sample run (the `e2e` golden tier against
      `Test/Sample.bam`), tracker cleanup, and now CI-observed-green (`refactor`
      pushed 2026-08-27; Lint run `33076602903` green against tip `d9327bf`).
- [x] `circos-removal/01` is done — all 8 boxes ticked; `src/ijump/circos.py`
      confirmed gone (no `circos*.py` anywhere in the tree).
- [x] `average-depth-zero-coverage/02` is done
- [x] `CHANGELOG.md` reads `## 2.0.0` (not `## Unreleased`) with all four breaking
      changes listed, including Circos removal — already true as of this check.
- [ ] `pyproject.toml`'s `version` is `2.0.0` — still reads `1.0.4`.
- [ ] `refactor` merged into `master` (merge commit or equivalent, per the
      maintainer's preferred git workflow)
- [ ] `master` tagged `v2.0.0`
- [ ] `ijump --help` and a manual real-sample run both work from a fresh
      checkout of the tag

## Comments
