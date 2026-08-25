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
