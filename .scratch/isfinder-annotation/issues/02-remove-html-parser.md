# 02 — Remove the ISFinder HTML parser

**What to build:** The HTML-scraping back-end ceases to exist. `ijump` no longer offers an
`isfinder-parse` subcommand, its help text lists the remaining three, and the README no
longer documents a workflow nobody can execute.

The ISFinder web site has been down for months, so its input cannot be produced. The
scraper also carried its own defects — a regex flag passed positionally where a count was
expected, so the flag never applied, plus the parsing problems recorded in
`review-followups 01`. Deleting it removes a maintenance surface rather than fixing one.

Independent of the rest of this series: it is a separate code path and touches no pipeline
output, so it does not need the goldens in place first.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] The HTML parser module and its test module are deleted
- [x] The CLI dispatch entry is removed; `ijump --help` lists exactly three subcommands
- [x] Invoking the removed subcommand name fails with the CLI's normal unknown-subcommand error
- [x] The README section describing the HTML workflow is removed
- [x] Full test suite green; no dangling import or reference to the deleted module anywhere in the repo

## Comments

- The deletion covers `src/ijump/isfinder_parser.py`, `tests/test_isfinder_parser.py`,
  the `cli.py` dispatch and help entries, the README's "2) Find them from ISFinder
  website" workflow (with its `Query=` header caveat), and the CLAUDE.md /
  `docs/agents/ast-grep.md` lines that named the module.
- `tests/test_cli_subcommands.py` replaces the deleted test module: it pins the
  rendered subcommand list to exactly the three survivors, and asserts the removed
  name now takes argparse's ordinary `invalid choice` path (exit 2).
- **`review-followups 01` closed `wontfix` as a consequence.** That ticket's
  remaining checkbox asked for `isfinder-parse` to be verified end to end under
  pandas 2.x, which this deletion makes unverifiable; its spec.md row and
  ticket 08's "blocked by 01" rationale were updated to match. Outside this
  ticket's stated scope, but it is the dangling tracker reference that criterion 5
  is about.
- Full suite: 69 passed, 4 skipped (e2e, inputs absent). `ruff check` /
  `ruff format --check` clean; `mypy` reports the same 6 pre-existing errors in
  `tests/test_clipped_read_search.py` as before the change.
