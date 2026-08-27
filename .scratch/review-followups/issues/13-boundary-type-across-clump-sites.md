# 13 — Give the recurring boundary clump a `Boundary` type

**What to build:** A `Boundary` dataclass (following `clipped_read_search.py`'s
existing `SearchResult` pattern) replaces the bare five-element lists and
positional unpacking at every one of its recurring call sites:
`_crtable_ungapped` (13 parameters), the frequency estimator (7), the
average-mode report (7), and the junction pairer (7). Callers read `.field`
instead of unpacking a list by position, and each function's parameter count
drops accordingly.

**Why:** review-followups/09 flagged this as the largest Data Clumps/Primitive
Obsession instance in the refactor. Scope was triaged (2026-08-26, see 09's
Comments) to all four remaining sites — the fifth, the Circos writer, is moot
because `.scratch/circos-removal/issues/01-remove-circos-module.md` deletes
that code entirely. Whether to also stop these functions mutating their
arguments in place was deliberately split out to review-followups/14, so this
ticket changes only the clump's shape, not the mutation semantics.

**Blocked by:** review-followups/12 (characterization tests must exist over
`_crtable_ungapped`'s path before its signature changes).

- [x] `Boundary` dataclass introduced, matching the five fields currently
      carried positionally
- [x] `_crtable_ungapped`'s signature threaded through `Boundary` where
      applicable, with its parameter count reduced accordingly
- [x] The frequency estimator, average-mode report, and junction pairer
      likewise take/return `Boundary` instead of raw positional lists —
      **scope corrected during implementation, see Comments**
- [x] No positional unpacking of boundary fields remains at any consuming site
      (of the sites that actually carry this shape)
- [x] In-place mutation semantics and the 3-tuple return shape are unchanged
      (that's ticket 14's job, not this one)
- [x] review-followups/12's characterization tests pass unchanged
- [x] Output byte-identical on a real sample before/after
- [x] `ruff`/`mypy` clean on `src/ijump/`

## Comments

**Scope correction (implementation time).** `find_pairs` (junction pairer),
`estimate_frequencies` (frequency estimator), and `report_average`
(average-mode report) do not actually take a `chrom+start+stop+edge+is_name`
group as parameters — `find_pairs` takes `(pos_l, pos_r, pos_l_count,
pos_r_count, chrom_len, max_is_dup_len, chrom)`; `estimate_frequencies` and
`report_average` take no chrom/start/stop/edge at all (those live in
dataframe columns, not parameters). None of the three have a `Boundary`-shaped
clump to fix. This ticket's "same shape recurs" language inherited 09's loose
framing without checking these signatures against the current code.

Confirmed with the maintainer (2026-08-26): apply `Boundary` only where the
five-field shape genuinely exists — `_crtable_ungapped`'s parameters, its
caller `clipped_read_search.search`'s `boundaries` list, and every other
consumer of that same list (`ISClipped.set_is_boundaries`,
`ISClipped._check_is_boundary_proximity`, the `backward_boundaries` list
comprehension in `ISClipped.run`). The other three functions are untouched by
this ticket. Ticket 14 (mutation removal) is scoped the same way: only
`_crtable_ungapped`, since that is the only function this ticket actually
changed.
