# 05 — Give the two duplicated definitions one home each

**What to build:** The empty summary-table shape and the minimum BLAST match
length are each defined once. Today each exists in two places, so changing one
and not the other produces output that disagrees with itself — silently, since
nothing compares them.

## Why

The module extractions left two things copied rather than shared.

**The summary-table initializer.** `ISClipped.sum_by_reg_tbl_init` and
`region_summary._sum_by_reg_tbl_init` have byte-identical bodies — same column
list, same `.extend` over the IS-element coordinate keys, same empty DataFrame.
**Both have live callers:** the `ISClipped` copy is used on the no-results path
to write the empty-output CSV header, and the `region_summary` copy is used
per-region while building the summary. So neither is dead code — this is a real
two-sources-of-truth problem. Add an IS element column to one and the empty
output file silently stops matching the populated one.

**The BLAST minimum.** `clipped_read_search.BLAST_MIN` is `10`, and
`ISClipped.blast_min` is also `10`. The extraction comment in
`clipped_read_search.py` records the intent honestly — *"Was ISClipped.blast_min;
never overridden by a caller"* — but the attribute was left behind on `ISClipped`
because it is **still live**: it is passed to `report_average` for the frequency
calculation. So the clipped-read search and the frequency report each read their
own copy of the same constant.

## Scope

- Make the empty-summary-table shape come from one definition that both the
  no-results path and the per-region path use. Watch the seam: the two current
  copies take their IS coordinates from different places (one from `self`, one
  from a parameter), so the shared version needs a parameter and the `ISClipped`
  side becomes a thin call.
- Give `BLAST_MIN` a single home that both `clipped_read_search` and the
  frequency-report path import, and delete the duplicate. `clipped_read_search`
  already owns the module-level constant with an explanatory comment, so that is
  the natural home unless importing it from there creates a cycle — check, and
  if it does, say where you put it instead and why.
- While there: `BLAST_MIN_IDENT` sits beside `BLAST_MIN` with the same
  "was an ISClipped attribute" comment. Confirm it has no surviving duplicate on
  `ISClipped` before assuming this ticket is done.

## Out of scope

- Changing either value. `10` and `75` stay `10` and `75`.
- Removing `ISClipped.blast_min`'s *use* — it legitimately feeds the frequency
  calculation. Only the second definition of the number goes away.
- The `min_match` parameter in the same file — that is ticket 04.

## Verification

- Output byte-identical on a real sample, both the populated path and the
  no-results path. The no-results path is the one at risk here and is the easy
  one to forget: `tests/test_no_results_paths.py` exercises it.
- Deliberately change the shared column list in a scratch edit and confirm
  **both** the empty output header and the populated summary change together —
  that is the actual property this ticket buys, and it is worth proving once
  before reverting the scratch edit.
- `tests/test_region_summary.py` and the no-results tests pass with no edits.
- `ruff`/`mypy` clean on `src/ijump/`; no import cycle introduced.

**Blocked by:** None — can start immediately. (Also edits `region_summary.py`,
which ticket 04 touches; different functions, so expect at most a trivial
conflict. Neither gates the other.)

**Status:** done

- [x] One definition of the empty summary-table shape, used by both callers.
- [x] One definition of the BLAST minimum, imported by both consumers.
- [x] `BLAST_MIN_IDENT` checked for a surviving duplicate.
- [x] Values unchanged; no import cycle.
- [x] Populated **and** no-results output byte-identical on a real sample.
- [x] Single-source property proven once by a scratch edit, then reverted.
- [x] Existing tests pass with no edits.

## Comments

Implemented. The shared table-init definition now lives in
`region_summary.py` as `sum_by_reg_tbl_init(is_coords)` (renamed from the
private `_sum_by_reg_tbl_init`, since it is now called across module
boundaries). `ISClipped.sum_by_reg_tbl_init` became a thin wrapper:
`return region_summary.sum_by_reg_tbl_init(self.is_coords)`.
`region_summary.summarize_by_region` was updated to call the renamed
function too.

`BLAST_MIN` kept its home in `clipped_read_search.py`
(`BLAST_MIN = 10`); `ISClipped.blast_min` now reads
`clipped_read_search.BLAST_MIN` instead of hardcoding `10`, preserving its
live use in `region_summary.report_average`'s frequency calculation while
eliminating the duplicate literal. No import cycle: `isclipped.py` already
imported both `clipped_read_search` and `region_summary` at module level,
and neither of those two modules imports `isclipped` back.

`BLAST_MIN_IDENT` checked via grep — no surviving duplicate on `ISClipped`
(only `self.blast_min` existed as a duplicate; `self.blast_min_ident` was
never present).

Scratch-edit proof performed as required: temporarily added a
`"SCRATCH_PROOF_COLUMN"` to the shared column list, confirmed both
`ISClipped.sum_by_reg_tbl_init()` (no-results path) and
`region_summary.sum_by_reg_tbl_init()` (populated path) picked up the new
column identically, then reverted.

Verification: `tests/test_region_summary.py` (6 tests) and the no-results
tests in `tests/test_no_results_paths.py` pass unchanged.
`tests/test_empty_run_outputs.py` (full CLI end-to-end, both average and
precise modes against `tests/fixtures/tiny.bam`) passes 4/4, proving the
no-results path's empty-header output is byte-correct post-refactor.
`ruff`/`mypy` clean on `src/ijump/`.

This ticket's `region_summary.py` change was merged alongside ticket 04's
(both touched the file; the merge additionally had to drop a leftover dead
`columns` local left behind at the merge site, consistent with ticket 04's
own removal of that same dead code — reconciled during the merge, not by
this ticket's implementing agent).
