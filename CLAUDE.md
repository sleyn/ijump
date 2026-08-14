# iJump

## Agent skills

### Issue tracker

Issues and specs live as markdown files under `.scratch/<feature-slug>/` in this
repo. See `docs/agents/issue-tracker.md`.

### Triage labels

The five canonical triage roles, each label string equal to its name. See
`docs/agents/triage-labels.md`.

### Domain docs

Single-context — one `CONTEXT.md` at the repo root plus `docs/adr/`. See
`docs/agents/domain.md`.

### Structural search

ast-grep for AST-based search and codemods, plus this repo's own rules in
`rules/` (run with `ast-grep scan .`). See `docs/agents/ast-grep.md`.
