# 01 — Precise mode's frequency formula was missing the blast_min correction

**What to build:** Precise mode corrects clipped-read counts for the same two
collection biases average mode corrects for, not just one.

## Why

Both modes' clipped-read evidence comes from the same BLAST-gated collection
(Step 2/3, `docs/algorithm.md`): a read's clipped segment must be ≥ 10bp
(`clipped_read_search.BLAST_MIN`) to ever be sent to BLAST at all, and the
read's matched (`M`) segment must be long enough for the aligner to place it
confidently. Average mode's `region_summary.report_average` corrects for
*both* biases — `(1 + blast_min/av_read_len)` adds back reads excluded by the
10bp BLAST gate, `(1 - min_match/av_read_len)` accounts for reads the aligner
couldn't place. `frequency_estimation.estimate_frequencies` (precise mode)
applied only the second correction. `blast_min` never appeared anywhere in
that file, even though precise mode's `N_cl`/`N_overlap` are built from the
identical 10bp-gated collection and suffer the identical undercount.

Found during a `/grilling` session (2026-08-25) auditing whether precise
mode's algorithm is scientifically correct, prompted by
`.scratch/average-depth-zero-coverage/issues/02-audit-frequency-correction-terms.md`
flagging average mode's correction terms as unreviewed — the same terms
recur in precise mode and were checked here.

## Scope

- Added a `blast_min_factor = 1 + blast_min/av_read_len` multiplicative
  factor to `N_clipped_l_correction`, `N_clipped_r_correction`,
  `N_overlap_l_correction`, `N_overlap_r_correction` in
  `frequency_estimation.py`, alongside the existing match-length factor —
  mirroring average mode's numerator correction exactly.
- `estimate_frequencies` gained a new `blast_min` parameter; the call site in
  `isclipped.py` passes `self.blast_min` (already
  `clipped_read_search.BLAST_MIN`).

## Out of scope

- `min_match`'s own robustness — see ticket 02.
- The wraparound pairing bug — see ticket 03, a separate function.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] `N_cl`/`N_overlap` corrections in precise mode apply the same two
      factors average mode applies, for the same reasons.
- [x] `tests/test_estimate_frequencies.py`'s golden values recomputed by
      running the function and re-pinned (not hand-derived).
- [x] `tests/goldens/e2e/precise/report/ijump_junction_pairs.txt` re-pinned;
      diff confined to the Frequency-derived columns, row count and pairing
      unchanged (44 rows before and after).
- [x] `ruff`/`ruff format`/`mypy` clean on `src/ijump/`.

## Comments

Full non-e2e suite (190 passed) and the full e2e tier (13 passed) both green
after the change. See tickets 02 and 03 in this same directory for the two
other findings from the same audit session, landed together.
