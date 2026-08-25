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

**Status:** ready-for-agent

- [ ] The HTML parser module and its test module are deleted
- [ ] The CLI dispatch entry is removed; `ijump --help` lists exactly three subcommands
- [ ] Invoking the removed subcommand name fails with the CLI's normal unknown-subcommand error
- [ ] The README section describing the HTML workflow is removed
- [ ] Full test suite green; no dangling import or reference to the deleted module anywhere in the repo
