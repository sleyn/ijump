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

**Status:** ready-for-agent

- [ ] README documents the IS-table columns and all three back-ends end to end
- [ ] No documentation references the removed HTML workflow
- [ ] Agent-facing architecture guidance describes the annotation stage and its back-ends
- [ ] Glossary terms for cluster, locus, fragment and origin-spanning are present and used consistently across the docs
- [ ] Both breaking changes are documented with the migration path for each
