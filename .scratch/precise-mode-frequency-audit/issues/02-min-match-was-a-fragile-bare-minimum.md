# 02 — min_match was a bare minimum over the whole run, not a robust statistic

**What to build:** The match-length correction coefficient is resistant to a
single outlier read, without changing what the statistic represents.

## Why

`min_match = min(match_lengths)` (in both `region_summary.report_average` and
`frequency_estimation.estimate_frequencies`) took the single smallest
matched-CIGAR-segment length pooled across every clipped read collected
anywhere in the run. One outlier — a spuriously short match from a marginal
alignment — sets the correction coefficient applied to every region's or
junction's reported frequency.

This is the min_match half of the open question in
`.scratch/average-depth-zero-coverage/issues/02-audit-frequency-correction-terms.md`
("do they use the right per-run aggregates ... or should they be scoped more
locally"). Investigated during the same `/grilling` session as ticket 01:
whether the fix should rescope the statistic locally (per-IS, per-region) or
keep it global. `min_match`/`av_read_len` model the *aligner's* own
minimum-placeable-match-length behaviour — a property of the aligner and the
read-length distribution, applied uniformly genome-wide by construction, not
a property of any individual IS element or region. Rescoping locally would
model something that isn't physically real, and would starve an already-noisy
per-read statistic down to the handful of reads supporting one junction. The
defect is the point estimator, not the scope.

## Scope

- Replaced `min(match_lengths)` with `np.percentile(match_lengths, 1)` in
  both `region_summary.py` and `frequency_estimation.py`, keeping the
  existing global, run-wide scope.

## Out of scope

- The `blast_min`/`av_read_len` term itself — untouched, and still an open
  question per `average-depth-zero-coverage/02`.
- Any change to how `match_lengths` is collected (still forward-search-only,
  still every clipped read in every IS boundary window).

**Blocked by:** None — can start immediately.

**Status:** done

- [x] `min_match` is a 1st-percentile estimate, not a bare minimum, in both
      modes.
- [x] `tests/test_region_summary.py` and `tests/test_estimate_frequencies.py`
      golden values recomputed by running the functions and re-pinned.
- [x] `tests/goldens/e2e/average/report/ijump_report_by_is_reg.txt` and
      `tests/goldens/e2e/precise/report/ijump_junction_pairs.txt` re-pinned;
      diffs confined to Frequency-derived columns, row counts unchanged.
- [x] `average-depth-zero-coverage/02` updated to record that its min_match
      half is resolved here, that precise mode shares the same computation
      and was affected by the same finding, and that the blast_min/av_read_len
      term remains open.

## Comments

Landed alongside tickets 01 and 03 in this same directory, from one audit
session. Full non-e2e suite (190 passed) and full e2e tier (13 passed) green.
