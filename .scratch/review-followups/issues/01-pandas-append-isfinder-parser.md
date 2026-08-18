# 01 — `isfinder-parse` crashes under the pandas version the recipe now resolves

**What to build:** `ijump isfinder-parse` runs to completion and writes its
output table in an environment built from the current `environment.yml` /
`meta.yaml`. Today it raises `AttributeError` partway through, because it calls
a pandas API that no longer exists in the version those files now allow.

## Why

`isfinder_parser.py` builds its result tables with three self-reassigning
`DataFrame.append` calls (grep `isels_ge = isels_ge.append`, `isels =
isels.append`, `isels_final = isels_final.append` — three sites, all in that one
module). `DataFrame.append` was deprecated in pandas 1.4 and **removed** in
pandas 2.0.

This was a known, accepted latent problem while the project pinned
`pandas=1.3.5`, where the method still worked. It is no longer latent:

- Packaging ticket 04 modernized the pin to `pandas<3` in both
  `environment.yml` and `meta.yaml`.
- Packaging ticket 02 wired this module up as the `isfinder-parse` subcommand of
  the `ijump` console script.

So a user creating the documented environment and running the documented
subcommand hits a hard crash. This is the highest-severity finding in the
review.

## Scope

- Rewrite the three `.append` call sites using `pd.concat`.
  `region_summary.py`'s `summarize_by_region` is the in-repo precedent — it
  already uses `pd.concat([...], sort=True)` for exactly this accumulate-rows
  pattern.
- **This is not a mechanical substitution.** The repo's own ast-grep rule says
  so and deliberately ships no autofix: `pd.concat` takes a *list*, its keyword
  arguments differ, and — the part that actually bites — appending to an
  **empty** DataFrame and concatenating with one differ in index and dtype
  handling. Check what each of the three accumulators is initialized to and
  confirm the resulting column order, index, and dtypes match what the old code
  produced.
- **Reconcile `rules/pandas-dataframe-append.yml`.** Its `note:` currently
  reads *"README.md pins pandas=1.3.5, where DataFrame.append still works but
  is deprecated. This is an upgrade blocker, not a live break."* Both halves are
  stale: README no longer pins a pandas version at all, and it **is** a live
  break. Update the note to match reality. Keep the rule itself and its
  `$X = $X.append($$$ARGS)` pattern — the rule's own note explains why the
  narrower self-reassignment shape is correct (the bare form returns 22 hits,
  21 of them ordinary list appends).

## Out of scope

- Re-litigating the `pandas<3` pin. It is correct; the code needs to catch up.
- Other pandas-2 incompatibilities. `ast-grep scan .` finds no others, so don't
  go hunting hypothetical ones.

## Verification

- `ast-grep scan .` reports zero `pandas-dataframe-append` findings.
- A test drives the IS-element parsing path on a small fixture and asserts the
  resulting table's shape, column order, and a couple of row values — enough to
  catch a `pd.concat` index/dtype regression, which is the actual risk here.
  Note `tests/` currently has no coverage of this module at all.
- The `isfinder-parse` subcommand runs end to end against a real ISFinder input
  under pandas 2.x and produces the same output as pandas 1.x did. If a
  pandas-1 comparison run isn't practical, say so in Comments rather than
  ticking this blind.
- `ruff check` / `ruff format --check` / `mypy` on `src/ijump/` still pass.

**Blocked by:** None — can start immediately.

**Status:** ready-for-agent

- [x] All three `DataFrame.append` sites rewritten with `pd.concat`.
- [x] Empty-accumulator index/dtype/column-order behaviour explicitly checked, not assumed.
- [x] `ast-grep scan .` clean for `pandas-dataframe-append`.
- [x] `rules/pandas-dataframe-append.yml`'s stale `note:` corrected.
- [x] A test covers the parsing path and would catch a `pd.concat` shape regression.
- [ ] `isfinder-parse` verified end to end under pandas 2.x. (Not practical to run here — see Comments.)
- [x] `ruff`/`mypy` still pass clean on `src/ijump/`.

## Comments

- **Empty-accumulator behaviour, checked not assumed.** Set up both a
  pandas==1.3.5 environment (via `uv run --no-project -p 3.10 --with
  "pandas==1.3.5" --with "setuptools<60" --with "numpy<1.22"`, the last
  pandas version where `.append()` still worked) and a pandas<3 environment
  (`uv run --no-project --with "pandas<3"`, matching this repo's current
  pin, resolved to 2.3.3) to compare directly. Ran the pre-rewrite (`git
  show HEAD:src/ijump/isfinder_parser.py`) and post-rewrite module against
  a hand-built ISFinder-HTML-shaped fixture (two IS elements on one contig,
  one on a second contig, one filtered out by the e-value cutoff) and
  diffed `ISTable_full.txt`/`ISTable_processing.txt` byte-for-byte. Also
  checked the zero-contig and zero-surviving-element edge cases the same
  way. All matched exactly, including a pre-existing quirk in the third
  site: `isels_final.append(isels, sort=False)` appends a frame indexed by
  `"name"` into a template where `"name"` is an ordinary column, so under
  old pandas the two never align — the `"name"` column ends up blank and a
  stray, all-blank `"check"` column survives (it was dropped from `isels`
  but is still in the template). Reproducing this exactly required seeding
  the `pd.concat` accumulator list with the original empty template (not
  just the collected per-contig frames), which is *not* what a naive
  "collect rows, concat once, skip if empty" rewrite would do — flagging
  this because it's the clearest instance of the ticket's "not a mechanical
  substitution" warning. Did not fix the quirk itself; out of scope per
  "Re-litigating the `pandas<3` pin... Other pandas-2 incompatibilities" and
  matching the spirit of ticket 06's DOTALL-bug precedent (flag, don't fix).
- **`isels_ge`/`isels` sites (the first two `.append` calls) rewritten
  differently from the third.** Both accumulate rows in a loop with no
  schema mismatch against their own template, so each collects rows into a
  plain list and does one `pd.concat` at the end, skipping the concat
  entirely when the list is empty (keeps the original empty, correctly
  columned template — also avoids pandas's "empty/all-NA entries"
  `FutureWarning` for those two sites, which a naive per-iteration
  `pd.concat([acc, new], ...)` substitution would trigger every time).
- **The third site (`isels_final`) does still emit that `FutureWarning`**,
  because it must include the empty template to reproduce the quirk above.
  This is expected and documented in a comment at the call site
  (`isfinder_parser.py`), not an oversight — code-review's Spec pass flagged
  that the warning-avoidance rationale on the first two sites could read as
  applying uniformly, so a clarifying comment was added directly above the
  third `pd.concat` call.
- **End-to-end verification against a real ISFinder BLAST HTML page under
  pandas 2.x was not practical in this environment**, so that checklist
  item is left unticked rather than ticked blind. Two blockers: (1) no real
  ISFinder output file exists anywhere in this repo or its `Example
  files`/`ijump_data` directories to use as input; (2) this sandbox cannot
  build the project's full dependency set at all (`pysam==0.15.4` /
  `pysamstats==1.1.2` fail to build from source here — `pkg_resources`
  build-dependency errors under both `uv sync` and manual
  `--no-build-isolation` workarounds, on Python 3.13/arm64 — matching the
  arm64/pysam wheel gap already documented in `README.md`'s Docker section).
  `isfinder_parser.py` itself only imports `pandas` (confirmed via `grep
  -n "^import\|^from" src/ijump/*.py`), so the rewrite was instead verified
  by (a) the byte-for-byte pandas-1-vs-pandas-2 comparison above using
  lightweight `uv run --no-project --with pandas<3`/`--with pandas==1.3.5`
  environments that don't need the rest of the stack, and (b) the new
  `tests/test_isfinder_parser.py`, run standalone the same way
  (`PYTHONPATH=src uv run --no-project --with "pandas<3" --with pytest
  pytest tests/test_isfinder_parser.py`). The full `pytest` suite could not
  be run in this sandbox for the same `pysam`/`pysamstats` build reason;
  ran the subset of test files with no `pysam`/`isclipped`/`sklearn`/`Bio`
  dependency instead (`tests/test_isfinder_parser.py`,
  `tests/test_estimate_frequencies.py`, `tests/test_find_pair.py`,
  `tests/test_region_summary.py`) — all 7 pass. `ruff check`, `ruff format
  --check`, and `mypy` were run scoped to `src/ijump/` and pass clean.
- **`docs/agents/ast-grep.md`** had a stale worked example (`Four hits
  today: isfinder_parser.py:55,89,115 and isclipped.py:1179... All four
  need hand-rewriting`) that this change makes doubly wrong: the three
  `isfinder_parser.py` sites are now rewritten, and the `isclipped.py:1179`
  reference was already stale on disk before this ticket (`grep -n append
  src/ijump/isclipped.py` shows only ordinary list `.append()` calls there
  now, resolved by an earlier isclipped-refactor ticket; `ast-grep scan .`
  is fully clean). Updated that bullet to reflect zero hits today. Outside
  this ticket's stated Scope, but it's a one-line factual-accuracy fix
  directly caused by this diff, so folded it in rather than filing a
  separate ticket for it.
