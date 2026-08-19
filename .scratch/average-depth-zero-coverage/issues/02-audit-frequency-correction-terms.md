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
