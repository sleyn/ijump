# 10 — Document the annotation subsystem

**What to build:** The cross-cutting documentation sweep, once all three back-ends and both
consuming modes are in place. Per-ticket documentation lands with each ticket; this one
covers what only makes sense described whole.

- The README's input section describes the IS table's columns and all three ways to produce
  one, and no longer describes the deleted HTML workflow.
- The architecture guidance for agents reflects the new module layout — one annotation
  stage with pluggable back-ends and a shared annotate-and-cluster core, replacing two
  standalone parser scripts.
- The domain glossary is complete and consistent, and the release notes state the two
  breaking changes plainly: the per-region report's column shape, and tables without a
  cluster column being rejected by precise mode.

**Blocked by:** 02, 05, 06, 07, 08, 09.

**Status:** done

- [x] README documents the IS-table columns and all three back-ends end to end
- [x] No documentation references the removed HTML workflow
- [x] Agent-facing architecture guidance describes the annotation stage and its back-ends
- [x] Glossary terms for cluster, locus, fragment and origin-spanning are present and used consistently across the docs
- [x] Both breaking changes are documented with the migration path for each

## Comments

**Release notes are a new `CHANGELOG.md`** — the repo had none. It leads with the two
breaking changes and the one-command remedy for each, then the additions, the three
frequency-reporting bugs this series fixed, and the removed HTML parser. Every number in it
was checked against the committed goldens rather than written from memory (the 0.98 case is
`N_clipped_r=284`, `N_unclipped_r=2`, `Frequency=0.9832…`). The README's opening points at
it.

**The HTML workflow criterion was already met** — the only surviving mentions of
`isfinder-parse` are in the test that asserts the removed subcommand is rejected, which is
the opposite of a stale reference. Verified across every `.md`, `.py`, `.toml`, `.yml` and
the Dockerfile.

**Three pieces of documentation were wrong rather than merely incomplete**, and are the real
content of this ticket:

- `docs/agents/domain.md` told agents "Python modules live flat at the repo root" and drew a
  tree of `isclipped.py`, `ijump.py`, `gff.py` at the root. The package moved to
  `src/ijump/` in the packaging work and this was never updated — so the one file that
  orients an agent before it explores was pointing at a layout that no longer exists.
- `docs/algorithm.md`'s Inputs table still described the IS table as four columns, and the
  write-up had no account of where that table comes from at all. It has an Annotation Stage
  section now: the three back-ends, the clustering rule with its thresholds and why coverage
  is measured on the shorter locus, cluster naming, and the origin-spanning flags. Every
  constant in it was checked against the code (95.0, 0.8, 1e-30, 0.75, 20).
- The README said legacy four-column tables are "still accepted", which is true of the
  reader and false of a run. It now says they are still *read* and points at the remedy.

**Smaller repairs.** The README's table of contents had not moved since before this series —
it listed neither the clusters section nor either new back-end nor the reports section — and
the Conda anchor was `<a name="'conda">`, with a stray apostrophe that had quietly broken
that link. A "Producing one" table now names all three back-ends in one place, which is what
"end to end" needed: they were each documented, but nowhere did the reader learn there were
three.

**Glossary.** `Fragment` was the missing term — used throughout the code and prose, defined
nowhere. Added, with the reason coverage is measured on the shorter locus. `Cluster`,
`Locus` and `Origin-spanning element` were already present and are used consistently;
`remnant` is left as accepted prose rather than rejected, since the code and README both use
it deliberately.

**Review follow-ups.** Both axes ran; between them they found one contradiction I introduced,
several inaccuracies, and three stale files this sweep had missed.

- **Two docs I wrote in the same commit disagreed.** `docs/algorithm.md` and `CLAUDE.md`
  called `isfinder-db-parse` "the only back-end that recovers family/group/pident from its
  input" while the README correctly said `isescan-convert` writes ISEScan's family. Only
  `group` and `pident` are ISFinder-only; corrected in both.
- **The CHANGELOG contradicted the README on which files carry the annotation stamp** — it
  implied `ijump_sum_by_reg.txt` does. It does not; the stamped pair is
  `ijump_report_by_is_reg.txt` and `ijump_junction_pairs.txt`.
- **A third breaking change was unlisted**: `isfinder-db-parse` now requires `-r/--ref`,
  because filling the cluster column means reading the loci out of the reference. Any
  existing invocation fails without it. Also added: the two threshold flags, and ticket 03's
  subject-id fix, which changes the emitted `is_name` for the eleven database entries whose
  own name contains an underscore.
- **Vocabulary drift the glossary should have caught.** The README's Clusters section said
  "two elements share a cluster", making a cluster a set of elements — inverting the
  definition. `CONTEXT.md` rejects "element" for a table row; those are loci. Fixed there,
  in the shared `--cluster-*` help text, and in `CONTEXT.md`'s own Cluster entry, which said
  "shorter element". The input list still called the IS table "File with mobile elements
  coordinates", a name the glossary explicitly rejects.
- **`cluster` meant two things in `docs/algorithm.md`** — the element grouping, and the
  agglomerative grouping of junction positions in precise mode. The second is "position
  cluster" now, with a note up front.
- **Three stale files this sweep had missed.** `docs/Average.md`'s worked example shows
  per-locus rows (`IS5_9`, `IS5_10`, …) that a current run collapses, and then walks the
  reader through summing them by hand — which is now what the tool does; it is framed as a
  pre-clustering run. `ijump_junctions.txt`'s *IS name* was documented as "mobile element
  name" in both mode docs when it is the locus the read matched. And
  `docs/agents/ast-grep.md`'s recipes still ran against root-relative `isclipped.py`, so
  they failed as written.
- **Two hardcoded counts in `docs/agents/ast-grep.md` had rotted** (38 vs 17, 22 vs 43).
  Rather than re-bake numbers that rot, the examples now state the *kind* of error and tell
  the reader to run both commands.

**Not done, needs filing.** The spec's out-of-scope note on the 75% overlap rule
("worth its own investigation, not a rider on this work") is described in `docs/algorithm.md`
as behaviour but has no ticket, unlike the union deferral which is durably recorded in the
README. `/to-tickets` is user-invocation only, so this needs one command from the maintainer.
