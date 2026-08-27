# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root
- **`docs/adr/`** — read ADRs that touch the area you're about to work in

If any of these files don't exist, **proceed silently**. Don't flag their absence; don't suggest creating them upfront. The `/domain-modeling` skill (reached via `/grill-with-docs` and `/improve-codebase-architecture`) creates them lazily when terms or decisions actually get resolved.

## File structure

This is a single-context repo. The package lives under `src/ijump/`; nothing importable sits
at the repo root.

```
/
├── CONTEXT.md
├── docs/
│   ├── adr/
│   │   ├── 0001-<slug>.md
│   │   └── 0002-<slug>.md
│   ├── agents/
│   └── algorithm.md
└── src/ijump/
    ├── cli.py               # dispatches every subcommand
    ├── ijump.py             # `ijump run` -- the detection pipeline
    ├── isclipped.py         # the pipeline's state and stages
    └── ...                  # see CLAUDE.md's Architecture section
```

`CLAUDE.md` is the authority on which module does what. Two shapes are worth knowing before
exploring, because both are easy to mistake for duplication:

- **The annotation stage** produces the IS table. Three back-ends read different inputs and
  all produce the same four locus columns, then hand off to one shared
  `is_annotation.annotate_and_cluster` that owns everything after. A rule about clusters or
  origin-spanning belongs in that core, never in a back-end.
- **The detection pipeline** consumes that table. `ISClipped` owns most of its state across
  methods — see `docs/agents/ast-grep.md`'s state-coupling matrix before extracting anything
  from it.

## Use the glossary's vocabulary

When your output names a domain concept (in an issue title, a refactor proposal, a hypothesis, a test name), use the term as defined in `CONTEXT.md`. Don't drift to synonyms the glossary explicitly avoids.

If the concept you need isn't in the glossary yet, that's a signal — either you're inventing language the project doesn't use (reconsider) or there's a real gap (note it for `/domain-modeling`).

## Flag ADR conflicts

If your output contradicts an existing ADR, surface it explicitly rather than silently overriding:

> _Contradicts ADR-0007 (event-sourced orders) — but worth reopening because…_
