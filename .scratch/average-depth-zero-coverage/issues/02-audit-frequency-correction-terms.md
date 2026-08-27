# 02 — Audit report_average's BLAST/match-length correction terms

Status: done
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

**Release decision (2026-08-26, release-2.0.0 planning).** The maintainer
chose to resolve the two remaining open questions now, before tagging 2.0.0,
rather than deferring — re-labelled `ready-for-agent` accordingly. This blocks
`.scratch/release-2.0.0/issues/01-cut-2-0-0-release.md`.

**Resolution (2026-08-26).**

*What bias `blast_min`/`av_read_len` corrects for.* A read whose true clip
point falls near either end of the read is systematically unobservable, for
two independent reasons:

- `clipped_read_search._write_cl_fasta` only writes a clipped segment to the
  BLAST query FASTA when it's `>= BLAST_MIN` (10bp,
  `clipped_read_search.py:25`). A read whose clip is shorter than that never
  becomes a BLAST hit, so it never enters `count`/`N_clipped_*` at all.
- Symmetrically, a read whose *matched* segment is very short is unlikely to
  have been placed reliably by the read aligner upstream of iJump in the
  first place — that's what `min_match` (the 1st percentile of observed
  match lengths, per `precise-mode-frequency-audit/02`) estimates.

Treating the true clip position as uniform across a read of length
`av_read_len`, only clip positions in the window
`[blast_min, av_read_len - min_match]` ever surface as an observed clipped
read. The observed count therefore undercounts the true count by
approximately `(av_read_len - blast_min - min_match) / av_read_len`, i.e. the
true count is `count * av_read_len / (av_read_len - blast_min - min_match)`.
The formula in the code applies two separate multiplicative factors instead
of that single ratio — `count * (1 + blast_min/av_read_len) /
(1 - min_match/av_read_len)`. Expanding both to first order in
`blast_min/av_read_len` and `min_match/av_read_len` gives the same
`1 + (blast_min + min_match)/av_read_len`, so the two-factor form is a
first-order approximation of the single-window correction, not an exact
match — accurate as long as `blast_min` and `min_match` stay small relative
to `av_read_len`, which holds for the short-read data iJump targets (10bp and
a low percentile of match length against reads on the order of 100-150bp).
No defect found; this is a deliberate, reasonable approximation, now spelled
out as a comment in `region_summary.report_average`
(`frequency_estimation.estimate_frequencies` already carried the equivalent
explanation from `precise-mode-frequency-audit/01`).

*`Depth = 0` / `NaN` interaction.* Already correctly handled and already
tested: `report_average` short-circuits to `Frequency = NaN` when
`Depth == 0` (`region_summary.py:108-119`) before the correction-term
division ever runs, so the `blast_min`/`min_match` factors can't produce
`inf`/`NaN`-from-division-by-zero surprises. Covered by
`test_report_average_frequency_is_nan_when_depth_is_zero` in
`tests/test_region_summary.py`, landed with ticket 01. Precise mode's
equivalent (`frequency_estimation.estimate_frequencies`) never divides by a
bare `Depth`; its `Frequency_l`/`Frequency_r` denominators carry a `+0.1`
pseudocount specifically to avoid a zero denominator, so no analogous NaN
case exists there.

Both open questions are resolved with no code-behavior change beyond an
explanatory comment; unblocks
`.scratch/release-2.0.0/issues/01-cut-2-0-0-release.md`.
