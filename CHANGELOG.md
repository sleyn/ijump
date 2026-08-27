# Changelog

## Unreleased

Dropped the `pysamstats` dependency; `average_depth`'s coverage calculation is now pure
`pysam` (a per-read CIGAR/span accumulator, no pileup), verified to reproduce `pysamstats`'
true (unrounded) coverage mean exactly across supplementary reads, internal
deletions/ref-skips, and zero-coverage windows.

**Fixed.** `average_depth` previously called `statistics.mean()` on a `numpy.int32` coverage
array, which truncated every fractional mean down to an integer. Coverage means (and, in
`--estimation_mode average`, the `Depth` column and everything derived from it) are now the
correct, unrounded float — if you compare a new run's `Depth`/`Frequency` values against one
from an earlier version, expect `Depth` to gain a fractional part it previously lost, and
`Frequency` (which divides by `Depth`) to shift slightly *downward* as a result, since the
old truncated, smaller `Depth` inflated it.

## 2.0.0

The IS table gained sequence-derived annotation, and both estimation modes now report by
**element** rather than by called locus. Four changes break existing setups; each has a
one-line remedy.

### Breaking: a run refuses a table with no `cluster` column

The `cluster` column groups the IS-table rows that are copies of one mobile element,
computed by aligning the called loci against each other rather than by comparing their
names. Both estimation modes group by it, so a table without it — a legacy four-column
table, or one with the column left blank — stops a run before it starts rather than being
guessed at.

**Migration.** Annotate the table you already have; its coordinates are preserved exactly:

```
ijump migrate-is-table -i <your IS table> -r <reference.fna> -d <ISFinder BLAST database> -o <outdir>
```

Or regenerate it with `ijump isfinder-db-parse`, which now writes the column itself.

### Breaking: `isfinder-db-parse` requires `-r/--ref`

Filling the `cluster` column means extracting each called locus from the reference and
comparing it with the others, so the parser now needs the genome the BLAST search was run
against. An existing `ijump isfinder-db-parse -b … -o …` invocation fails until `-r` is
added; nothing else about it changed.

### Breaking: the per-region report names elements, not loci

`ijump_report_by_is_reg.txt` and `ijump_sum_by_reg.txt` carry one entry per cluster where
they used to carry one per IS-table row. On the test genome the three entries `IS17_1`,
`IS17_2` and `ISAba12_1` become the single entry `ISAba12` — they are one copy and two of
its own fragments, and splitting one insertion's evidence three ways understated every
frequency and could hide an insertion under the reporting cutoff.

**This changes numbers, not just names.** Where two loci of one element hit the same region,
their counts now sum.

**Migration.** Any script selecting rows by IS name (`IS17_1`) needs the cluster name
(`ISAba12`) instead. `ijump combine-results` handles it for you. The names are in the IS
table's `cluster` column, so a lookup from old name to new is one read of that file.

`ijump_report_by_is_reg.txt` — and precise mode's `ijump_junction_pairs.txt`, the two files
`combine-results` merges — also gained a leading `# ijump-is-table: <digest>` line naming
the annotation the run used. `ijump_sum_by_reg.txt` is not merged and is not stamped. `combine-results` reads it and **refuses to merge
samples annotated against different IS tables**, which would otherwise line up cluster names
that mean different elements. Reports written before this change carry no such line and are
refused; rerun those samples.

### Breaking: Circos output is gone

`-c/--circos` and `circos.py` are removed. Circos rendering had no role in detection itself —
it only turned finished detection results into Circos input files — and saw little practical
use. An invocation passing `-c/--circos` now fails with an unrecognized-argument error; drop
the flag and, if you relied on the generated Circos files, render them from
`ijump_report_by_is_reg.txt`/`ijump_junction_pairs.txt` with your own tooling instead.

### Added

- `ijump migrate-is-table` — annotate an existing IS table instead of regenerating it.
- `ijump isescan-convert` — build an IS table from ISEScan results. iJump reads ISEScan
  output and never runs ISEScan. The two annotations are complementary: ISEScan finds
  elements with no ISFinder hit at all, while ISFinder finds fragments too short for
  ISEScan's terminal-repeat-plus-ORF test.
- `family`, `group`, `pident`, `wraps_origin` and `element_id` columns in the IS table.
  `wraps_origin` and `element_id` mark a copy the assembler cut at a contig boundary.
- `--cluster-identity` and `--cluster-coverage` on every IS-table back-end, for the
  clustering thresholds (default 95% identity over 80% of the shorter locus). Both reject a
  value off their own scale rather than answering wrongly: `--cluster-coverage 80` is an
  error, not a request for 80× coverage.

### Fixed

- A cluster whose junctions on a contig were all right-handed had them written into the
  left-hand columns, so frequency estimation read its evidence from the wrong side and
  reported `0.0`. On the test sample one such insertion had 284 clipped reads against 2
  unclipped — a frequency of 0.98 reported as absent.
- Position 0, the first base of a contig, was indistinguishable from "this junction has no
  partner", so a junction there was dropped or silently zeroed.
- Eleven ISFinder database entries carry an underscore inside the element's own name
  (`ISBj2_B_IS5_IS5`). Splitting the subject id from the left truncated those to `ISBj2`,
  so the IS table named the locus after the wrong element. The split runs from the right
  now, which changes the emitted `is_name` for any genome carrying one of the eleven.
- `combine-results`' lab format derived its element name by stripping a numeric suffix,
  which truncates a cluster name instead: `ISAba1`, `ISAba11` and `ISAba12` all collapsed to
  `ISAba`, summing unrelated elements into one row.

### Removed

- The ISFinder HTML parser, whose upstream web site went down. Use the BLAST workflow
  (`ijump isfinder-db-parse`) instead.

## 1.0.4

### Fixed

- `--estimation_mode precise` was compared against the misspelled `'presice'`, so the `IS
  pos` column of `ijump_junctions.txt` in precise mode was left 0-based while `Position` on
  the same row was already 1-based. `IS pos` is now converted to 1-based too, so if you
  compare a new run's `ijump_junctions.txt` against one produced by an earlier version, `IS
  pos` will be shifted by 1.
