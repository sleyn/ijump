# 14 — Remove in-place mutation from the Boundary-touched functions

**What to build:** `_crtable_ungapped` and the other functions touched by
review-followups/13 stop mutating their `Boundary`/list arguments in place and
instead return everything they compute explicitly. Callers no longer need to
know that a value they pass in will be changed out from under them.

**Why:** review-followups/09's triage (2026-08-26) confirmed the mutation
cleanup is wanted, but explicitly as a separate, later change from the
`Boundary`-type introduction — doing both signature shape and mutation
semantics in one pass was called out as the risky combination to avoid.

**Blocked by:** review-followups/13 (the `Boundary` type must exist first;
this ticket changes how it's passed, not what it looks like).

- [x] None of the touched functions mutate a parameter in place; all outputs
      are returned explicitly
- [x] Call sites updated to use the returned values instead of relying on
      mutation
- [x] review-followups/12's characterization tests (extended if needed to
      assert non-mutation) pass
- [x] Output byte-identical on a real sample before/after
- [x] `ruff`/`mypy` clean on `src/ijump/`

## Comments

Scoped to `_crtable_ungapped` only, per 13's scope correction (see that
ticket's Comments) — it is the only function 13 actually touched.
`_crtable_ungapped` now builds `clipped_reads`/`cl_read_cov_overlap`/
`match_lengths` as locals scoped to one `Boundary` and returns all six values
(`index, read_lengths, n_reads_analyzed, clipped_reads, cl_read_cov_overlap,
match_lengths`); its only caller, `clipped_read_search.search`, merges each
call's return into its own running totals (`dict.update` for clipped reads,
whose index ranges never overlap between boundaries; a summing merge for
`cl_read_cov_overlap`, since the same reference position can recur across
boundaries and the original mutated a single shared dict with `+=`).
`_cl_read_cov_overlap` (the private helper `_crtable_ungapped` calls) is
unchanged — it still mutates the dict it's handed, but that dict is now
`_crtable_ungapped`'s own local, never a reference any caller holds onto.
