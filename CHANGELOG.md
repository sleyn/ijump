# Changelog

## Unreleased

The IS table gained sequence-derived annotation, and both estimation modes now report by
**element** rather than by called locus. Two changes break existing setups; both have a
one-command remedy.

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

The two files `combine-results` merges also gained a leading `# ijump-is-table: <digest>`
line naming the annotation the run used. `combine-results` reads it and **refuses to merge
samples annotated against different IS tables**, which would otherwise line up cluster names
that mean different elements. Reports written before this change carry no such line and are
refused; rerun those samples.

### Added

- `ijump migrate-is-table` — annotate an existing IS table instead of regenerating it.
- `ijump isescan-convert` — build an IS table from ISEScan results. iJump reads ISEScan
  output and never runs ISEScan. The two annotations are complementary: ISEScan finds
  elements with no ISFinder hit at all, while ISFinder finds fragments too short for
  ISEScan's terminal-repeat-plus-ORF test.
- `family`, `group`, `pident`, `wraps_origin` and `element_id` columns in the IS table.
  `wraps_origin` and `element_id` mark a copy the assembler cut at a contig boundary.

### Fixed

- A cluster whose junctions on a contig were all right-handed had them written into the
  left-hand columns, so frequency estimation read its evidence from the wrong side and
  reported `0.0`. On the test sample one such insertion had 284 clipped reads against 2
  unclipped — a frequency of 0.98 reported as absent.
- Position 0, the first base of a contig, was indistinguishable from "this junction has no
  partner", so a junction there was dropped or silently zeroed.
- `combine-results`' lab format derived its element name by stripping a numeric suffix,
  which truncates a cluster name instead: `ISAba1`, `ISAba11` and `ISAba12` all collapsed to
  `ISAba`, summing unrelated elements into one row.

### Removed

- The ISFinder HTML parser, whose upstream web site went down. Use the BLAST workflow
  (`ijump isfinder-db-parse`) instead.
