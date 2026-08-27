# 06 — Naming sweep: two typos, one unfinished rename, one glossary breach

**What to build:** Someone grepping this codebase for a concept finds it under
one spelling. Four small naming defects landed during the extraction work; each
is individually trivial, and together they mean grep no longer finds what it
should.

## Why

All four came in with the refactor batch, and none is caught by tooling —
ruff and mypy are perfectly happy with a misspelled name used consistently.

1. **`_resore_orig_counts`** in `frequency_estimation.py` — missing a `t` in
   "restore". A new module-level function, so the typo is now the public-ish
   name of that step, referenced at two call sites and in a comment.
2. **`_cirocs_colors`** in `circos.py` — "circos" transposed, in a module
   literally named `circos`. Referenced at roughly six sites.
3. **An unfinished rename.** isclipped-refactor ticket 12 renamed
   `ISClipped._av_depth` to `average_depth` specifically to close this seam, but
   the circos extraction kept the old spelling for its parameter: `write_files`
   still takes `av_depth`. Note this parameter is a **passed-in callable**, not a
   scalar depth value, which is worth reflecting in whatever name it gets.
   `circos.py`'s own header comment already flags this parameter as a late
   addition outside ticket 11's enumerated signature.
4. **A glossary breach.** `CONTEXT.md`'s Language section defines **Region** —
   *"An annotated genetic element — a gene or an intergenic interval — that
   junctions get grouped into for average-mode reporting"* — and its `_Avoid_`
   list names **GE** explicitly. Three comments in `isclipped.py` use
   "genetic element (GE)" and "each GE". Two of them sit directly above the
   average-mode report calls, i.e. exactly where the documented term matters
   most.

## Scope

- Fix both typos and all their references.
- Finish ticket 12's rename through the circos parameter and its uses inside
  that module. Since it is a callable, consider whether the name should say so.
- Replace "GE"/"genetic element (GE)" with the glossary term in comments.
  **Re-grep for other occurrences** across `src/ijump/` rather than fixing only
  the three the review found — including the `_Avoid_` list's other entries
  ("locus", and "gene" where a region may be intergenic).
- Rename references in comments and docstrings too, not just code identifiers.
  A comment that still says `av_depth` is the same grep failure.

## Out of scope

- Renaming anything user-facing: output column headers, output filenames, CLI
  flags, config keys. The glossary governs *our* vocabulary, not the file format
  users already parse. If a documented term collides with an output header, leave
  the header alone and note it in Comments.
- Restructuring `write_files`' signature — it has too many parameters, but that
  is ticket 09's question. This ticket renames one of them and stops.
- `docs/` prose. If the glossary is breached there too, note it; don't expand.

## Verification

- `grep -rn` for each old spelling across `src/` returns nothing.
- Output byte-identical on a real sample. These are pure renames — any output
  change means something was misidentified as a rename when it wasn't.
- Full test suite passes with no test edits. If a test referenced an old name,
  updating it is fine, but note it — it means the name was more public than
  assumed.
- `ruff`/`mypy` clean on `src/ijump/`.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] `_resore_orig_counts` corrected, all references updated.
- [x] `_cirocs_colors` corrected, all references updated.
- [x] The circos depth-callable parameter renamed, finishing ticket 12.
- [x] "GE" replaced with the glossary term; `src/ijump/` re-grepped for other `_Avoid_`-list breaches.
- [x] Comments and docstrings updated alongside identifiers.
- [x] No user-facing name changed.
- [x] Output byte-identical on a real sample.

## Comments

Implemented. `_resore_orig_counts` → `_restore_orig_counts` in
`frequency_estimation.py` (definition, 2 call sites, 1 comment).
`_cirocs_colors` → `_circos_colors` in `circos.py` (definition + 4 usage
sites).

`write_files`'s `av_depth` parameter renamed to **`average_depth_fn`**
rather than bare `average_depth`, because `region_summary.report_average`
already has a parameter literally named `average_depth` for the same kind
of callable — reusing that name in `circos.py` would relocate rather than
close the scalar/callable ambiguity ticket 12 set out to fix. Header
comment and both internal call sites updated; the external caller
(`ijump.py`) passes positionally so needed no change.

Replaced "genetic element (GE)"/"each GE" with "Region" at all 3
`isclipped.py` comment sites. Re-grepped `src/ijump/` against every
`CONTEXT.md` `_Avoid_` entry (GE, locus, gene, empty result/no data/
failure/error, approximate/fast/exact/detailed mode) — no other breaches
found; remaining `gene`/`Locus tag` hits are genuine GFF-annotation fields
or user-facing output column headers, correctly left alone (out of scope
per the ticket).

No test edits — `tests/test_circos.py` still says `av_depth` in its
docstring and `fake_av_depth` fixture name, but the call there is
positional so it needed no change; this is a residual outside the
ticket's `src/`-scoped verification.

Verification: `grep -rn` for all three old spellings across `src/` returns
nothing. `ruff`/`mypy` clean on `src/ijump/`. Output verified
byte-identical on `tests/fixtures/tiny.*` before/after. Full suite: 35
passed / 11 failed in the implementing sandbox, confirmed via `git stash`
that the same 11 tests fail identically with and without this diff — a
pre-existing, documented `pysamstats` environment gap (ticket 16),
unrelated to this change.
