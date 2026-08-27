# 13 — Extract average mode's report generation into a standalone `region_summary` module

Status: ready-for-agent
Blocked by: —

## Why

`summary_junctions_by_region` (`isclipped.py:637-664`) and `report_average`
(`isclipped.py:676-707`) are average mode's whole algorithmic core — the counterpart to
precise mode's frequency-estimation stage, which ticket 09 already extracted into
`frequency_estimation.py`. Average mode has no equivalent module: its report logic sits
directly on `ISClipped`, reaching into `self.junctions`, `self.is_coords`, `self.gff.ann_pos`,
`self.match_lengths`, `self.read_lengths`, `self.n_reads_analyzed`, `self.blast_min`,
`self.min_match`, and `self.average_depth`.

It is also the only remaining corner of `ISClipped` with no direct test coverage —
`tests/test_isclipped.py` only exercises `average_depth`. `summary_junctions_by_region`
(`isclipped.py:659`) calls `self.sum_by_region = self.sum_by_region.append(temp, sort=True)`,
flagged by this repo's own `rules/pandas-dataframe-append.yml`: not a live break (pandas is
pinned to 1.3.5 per `README.md:71`, where `.append` still works, just deprecated), but removed
in pandas 2.0 and untested either way — nothing currently exercises the branch where a second
IS element lands in the same annotated region.

Architecture review candidate 01 (`/improve-codebase-architecture`, 2026-08-13), grilled to a
shared understanding the same session. Vocabulary: `/codebase-design`.

Decisions from grilling, all load-bearing for the scope below:

1. **`average_depth` stays on `ISClipped`.** Ticket 12 settled this the same day: it's a lazy
   callback shared between `report_average` and `circos.write_files`, and its `@lru_cache`
   dedupes `pysamstats` calls across both callers. Moving it here would re-litigate ticket 12.
   It is passed into `report_average` as an injected callable, same pattern `circos.write_files`
   already uses.
2. **Two public functions, not one.** `ISClipped.run()` writes `sum_by_region` to CSV
   (`isclipped.py:556`) *before* computing `report_table` from it (`isclipped.py:559-560`) —
   two independently-written outputs, not one bundled result. No dataclass wrapper.
3. **`report_average`'s `min_match`/`av_read_len` mutation of `self` is dropped.** Confirmed by
   grep: nothing outside `report_average` reads `self.min_match` or `self.av_read_len` after
   it returns (both are read at `isclipped.py:701-702`, inside the same method that sets them).
   They become local variables in the extracted function.
4. **Fix the `.append()` call in this ticket**, not a separate one — it's the exact method
   being relocated, the rewrite is one line (`pd.concat`), and it's the last item this repo's
   pandas-append rule flags in `isclipped.py`.

## Scope

### 1. Create the module

Create `region_summary.py` at the repo root (alongside `junction_pairing.py`,
`frequency_estimation.py`, `clipped_read_search.py`) with two public functions:

```python
def summarize_by_region(junctions, is_coords, gff_ann_pos) -> pd.DataFrame:
    ...

def report_average(sum_by_region, match_lengths, read_lengths, n_reads_analyzed,
                    blast_min, min_match, average_depth) -> pd.DataFrame:
    ...
```

- `summarize_by_region` moves the body of `summary_junctions_by_region`
  (`isclipped.py:637-664`), reading `junctions`/`is_coords`/`gff_ann_pos` as parameters instead
  of `self.junctions`/`self.is_coords`/`self.gff.ann_pos`. Replace
  `self.sum_by_region = self.sum_by_region.append(temp, sort=True)` (`isclipped.py:659`) with
  `pd.concat([sum_by_region, temp], sort=True)` (or equivalent) — this is the one intentional
  behavior change in this ticket, see Verification.
- `report_average` moves the body of `report_average` (`isclipped.py:676-707`), taking
  `average_depth` as a callable parameter (called the same way `self.average_depth(chrom,
  start, stop)` is today) instead of reaching into `self`. `min_match` and `av_read_len`
  become local variables, computed exactly as today (`min(match_lengths)`,
  `read_lengths / n_reads_analyzed`) — not assigned back onto anything (see "Why," decision 3).

Preserve each body otherwise verbatim — relocation plus the one named fix, not a rewrite.

### 2. Update the call site

`ISClipped.run()`'s average branch (`isclipped.py:545-556`, `558-560`):

```python
self.sum_by_region = region_summary.summarize_by_region(
    self.junctions, self.is_coords, self.gff.ann_pos
)
self.sum_by_region.to_csv(os.path.join(outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)

self.report_table = region_summary.report_average(
    self.sum_by_region, self.match_lengths, self.read_lengths, self.n_reads_analyzed,
    self.blast_min, self.min_match, self.average_depth,
)
self.report_table.to_csv(os.path.join(outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)
```

Add `import region_summary` to `isclipped.py`. Remove `summary_junctions_by_region` and
`report_average` from `ISClipped` once moved.

## Out of scope

- `average_depth` (`isclipped.py:667-673`) — stays on `ISClipped` (decision 1).
- Deduplicating `min_match`/`av_read_len` against `frequency_estimation.estimate_frequencies`'s
  independent computation of the same formulas (noted, not touched, by ticket 09). Still not
  touched here.
- `ISClipped.run()`'s broader shape (splitting average/precise orchestration into separate
  methods) — a separate architecture review candidate, not this ticket.

## Verification

**Characterization test — single-hit-per-region path.** Write
`tests/test_region_summary.py` against *today's* `summary_junctions_by_region`/
`report_average` for the case exercised in practice (each annotated region receives at most
one junction per IS element — the `else` branch at `isclipped.py:660-662` is never hit).
Hand-construct a small `junctions`/`is_coords`/`gff_ann_pos` fixture, pin today's output as a
golden value before the move, confirm the moved functions reproduce it after.

**Fresh spec-based test — multi-hit-per-region path.** No trustworthy current behavior exists
to pin for the case `.append`'s replacement changes (two-or-more junctions for the same IS
element landing in the same region) — nothing exercises it today. Write a new test asserting
`summarize_by_region` correctly accumulates counts across multiple hits in one region, against
the `pd.concat`-based implementation directly. State in the test docstring that this path has
no historical baseline and the assertion is spec-based, not characterization.

**Manual real-sample diff — single-hit path only.** Run the full pipeline in average mode on a
real sample before and after this change, diff `ijump_sum_by_reg.txt` and
`ijump_report_by_is_reg.txt` byte-for-byte. Must match exactly (real samples seen so far don't
hit the multi-hit branch, consistent with `.append` never having been exercised there).

Command to use, both before and after:

```
python ijump.py \
  -a ./Test/Sample.bam \
  -r ./Test/A_baumannii_assembly.fna \
  -g ./Test/Acinetobacter_baumannii_ATCC17978.gff \
  -i ./Test/ISTable_processing.txt \
  -o ./Test/Junctions_average/
```

Copy `./Test/Junctions_average/ijump_sum_by_reg.txt` and `ijump_report_by_is_reg.txt` aside
before making the change (e.g. to `ijump_sum_by_reg.txt.before`/
`ijump_report_by_is_reg.txt.before`), rerun the same command after, and diff.

## Done when

- `region_summary.py` exists with public `summarize_by_region(...)` and `report_average(...)`
  functions; `ISClipped` no longer defines either method's body.
- `isclipped.py:659`'s `.append()` call no longer exists anywhere in the codebase.
- `ISClipped.run()` calls `region_summary.summarize_by_region` and `region_summary.report_average`
  and assigns their results to `self.sum_by_region`/`self.report_table`.
- `tests/test_region_summary.py` exists with both the single-hit characterization test and the
  multi-hit spec-based test, and both pass.
- Manual real-sample before/after diff of `ijump_sum_by_reg.txt` and `ijump_report_by_is_reg.txt`
  is clean.
- `pytest` passes from a clean clone.

## Comments
