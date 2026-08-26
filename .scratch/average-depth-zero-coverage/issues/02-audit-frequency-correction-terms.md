# 02 — Audit report_average's BLAST/match-length correction terms

Status: needs-triage
Blocked by: —

## Open question

`region_summary.report_average`'s `Frequency` column (region_summary.py:85-92)
applies two correction terms on top of the raw junction-count ratio:

```python
report_table["Frequency"] = report_table.apply(
    lambda x: round(
        (x["count"] / 2 * (1 + blast_min / av_read_len))
        / (x["Depth"] * (1 - min_match / av_read_len)),
        4,
    ),
    axis=1,
)
```

where `min_match = min(match_lengths)` and `av_read_len = read_lengths / n_reads_analyzed`.

During a grilling session on `average_depth` and `report_average`
(2026-08-17), the `count / 2` term was walked through and confirmed correct
(it corrects for one genomic insertion event producing junction-supporting
reads at both the IS element's left and right boundaries). The
`blast_min`/`min_match`/`av_read_len` correction terms were explicitly *not*
reviewed — the session stopped before reaching them, flagged as a separate
sampling-bias question deserving its own focused investigation.

## Scope

No fix is specified here. This ticket exists to track that the correction
terms have not yet been validated:

- What sampling bias are `(1 + blast_min / av_read_len)` and
  `(1 - min_match / av_read_len)` meant to correct for?
- Do they use the right per-run aggregates (`min` of all match lengths,
  overall average read length) or should they be scoped more locally (e.g.
  per-region, per-IS)?
- Are they still correct after ticket 01's `Depth = 0` / `Frequency = NaN`
  change?

Needs triage/scoping before it's ready for an agent to implement.

## Comments

**The `min_match`/`av_read_len` scoping half is resolved —
`.scratch/precise-mode-frequency-audit/issues/02-min-match-was-a-fragile-bare-minimum.md`
(2026-08-25).** `min_match`/`av_read_len` model the *aligner's* own
minimum-placeable-match-length behaviour — a property of the aligner and the
read-length distribution, applied uniformly genome-wide by construction, not
a property of any individual IS element or region — so the answer to "should
they be scoped more locally" is no: local scoping would model something that
isn't physically real, and would starve an already-noisy per-read statistic
down to the handful of reads supporting one junction or region. The real
defect was the point estimator: `min()` let a single outlier read set the
correction applied to every region in the run. Replaced with the 1st
percentile of `match_lengths`, keeping the global scope. Landed in both
`region_summary.report_average` (this function) and
`frequency_estimation.estimate_frequencies` — precise mode's frequency
formula does the identical `min_match`/`av_read_len` computation and was
affected by the same finding (`precise-mode-frequency-audit/01` also found
precise mode missing the `blast_min` term entirely, a separate defect).

**Still open**: the `blast_min`/`av_read_len` sampling-bias question itself
(what bias each term is meant to correct for, beyond "adds back reads
excluded by the BLAST/aligner gates") and the third bullet (`Depth = 0` /
`NaN` interaction) — neither was investigated in the 2026-08-25 session.
