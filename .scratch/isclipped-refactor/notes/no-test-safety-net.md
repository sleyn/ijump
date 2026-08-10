# Note: there is no test safety net for the ISClipped refactor

Status: note (superseded by ticket 06)
Found: 2026-08-05, during setup, before `/improve-codebase-architecture` ran

## The finding

`Test/test_find_pair.py` is the only test touching `ISClipped`'s hardest function,
and it verifies nothing.

Two independent defects:

1. **It tests a copy, not the code.** The file pastes its own definition of
   `_find_pair` at line 5 rather than importing from `isclipped.py:667`. The
   shipped implementation is never loaded. The copy can drift from the real one
   silently — and given the two are now four years apart in provenance, it likely
   already has.

2. **It has no assertions.** Line 173 calls the function and binds the result to
   `x`; line 174 prints `'OK'` unconditionally. The script exits 0 whatever the
   function returns. Verified by running it after moving the stale modules out —
   still prints `OK`.

`Test/iJump_test/test_pysam_pileups.py` (36 lines) is the only other test file.

## Why it matters for the refactor

`_find_pair` is ~155 lines (`isclipped.py:667-821`) and is the junction-pairing
core — the part of `ISClipped` most likely to change shape under refactoring and
hardest to reason about by eye. Right now a refactor of it cannot be validated by
anything.

More broadly: `isclipped.py` is 1333 lines across 35 methods on one class, with a
dozen DataFrames as shared instance state, and effectively zero executable
verification. Any extraction is currently unfalsifiable.

## What to do with it

Not a delete, and not a quick fix — this is the seed of the characterization
harness that should be **ticket #1** of the refactor, ahead of any extraction.

Usable material already present in the file:

- A real input fixture at lines 165-171 — `pos_l`, `pos_r`, `pos_l_count`,
  `pos_r_count`, `chrom_len = 3909467`, `max_is_dup_len = 20`, and a real contig
  name. Captured from actual data, so it exercises real coordinate ranges.

Shape of the work:

- Import the real `_find_pair` from `isclipped.py`; delete the pasted copy.
- Pin current output as a golden value and assert on it. Characterization, not
  specification — the goal is "did I change behaviour", not "is this correct".
- Extend beyond `_find_pair` to an end-to-end snapshot over the sample data in
  `Test/` (`Sample.bam`, the A. baumannii references) so the full pipeline has a
  baseline too.

## Constraint worth carrying

`Test/` is gitignored (`.gitignore:4`) and holds ~8.7 GB of BAMs and BLAST
databases. Test *code* needs to move somewhere tracked; the fixture *data* cannot
be committed as-is. Deciding where tests live and how they reach their data is
itself a design question — likely worth raising during
`/improve-codebase-architecture` or the grilling that follows.

## Related cleanup already done

`Test/isclipped.py` (`class isclipped`, pre-`direction` API) and `Test/gff.py`
(still using `is '+'` identity comparison) were stale ancestors, imported by
nothing. Moved to `~/ijump-attic/` rather than deleted, since `Test/` is
gitignored and they were never in version control. Delete for real once the
refactor is underway and their absence has clearly caused no problem.
