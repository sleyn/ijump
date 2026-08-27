# 07 — Average mode reports one column per cluster

**What to build:** The per-region report emits one column per cluster instead of one per IS
name. On the reference genome that turns three columns — `IS17_1`, `IS17_2`, `ISAba12_1`,
which are one element and two of its own fragments — into a single column for the element
that actually jumped.

**This changes the shape of a documented collectable output.** Anyone with downstream
scripts reading the per-region report by column name is broken by it. That is the point of
the change and it should be called out in the release notes, not smoothed over.

Includes the open risk recorded in the spec: the multi-sample merge joins per-sample
reports on IS identity, which is only sound when every sample in a merge was annotated
against one reference and one IS table. Cluster names are *derived* rather than fixed, so
merging samples annotated against different references would silently misalign them —
strictly worse than today. The merge needs to detect that and refuse.

**Blocked by:** 04 — Similarity clusters computed and written to the table.

**Status:** done

- [x] Per-region reporting emits one column per cluster, named for the cluster
- [x] The per-region report golden is re-pinned, with the three `IS17`/`ISAba12` columns collapsing to one
- [x] The multi-sample merge detects reports built from different IS tables and refuses with a clear message
- [x] Merging reports that do share an IS table behaves as before
- [x] The output-shape change is documented as breaking

## Comments

**The collapse is bigger than the ticket's example, and it changes numbers.** `ISAba11_1/2/3`
collapse to `ISAba11` for the same reason `IS17`/`ISAba12` do, and where two loci of one
element hit the same region their counts now sum into one row: the first report row goes
from frequency 0.0017 to 0.0034. That means the diff is not a relabelling.

The clearest case is region `AUO97b_01699` at 1658711. It was two rows, `IS17_1` and
`IS17_2`, at 0.0043 each; it is now one `ISAba12` row at 0.0085. Circos draws at a cutoff of
0.005, so neither half cleared it and the insertion was **absent from the diagram
entirely** — it appears in `links.txt`, `text.txt` and `histogram.txt` for the first time
here. (It was always present in the per-region report itself, at half strength in each of
two rows. An earlier draft of this comment said it was below the *report's* cutoff; that was
wrong, and the review caught it.)

**Circos had to follow, or the run crashes.** `write_files` looked up `is_coords[is_name]`
with a name taken from the report, which is now a cluster and not a row of the IS table, and
its histogram indexed the region summary by locus name. It takes the cluster lookup now:
labels stay per-locus (that is where the elements physically are), while colour and the
histogram's columns follow the cluster, and a link is drawn from *every* locus of the
cluster — any copy could be the one that jumped, which is precisely why they share a
cluster. Not in the ticket, not optional.

**Average mode now requires the cluster column**, which ticket 06 had deliberately left
alone ("that is isfinder-annotation 07's call to make"). The test asserting average mode
accepts a legacy table is flipped to assert it refuses.

**How the merge detects mixed annotations.** Nothing in a report's columns answers "which IS
table produced this", and comparing the sets of IS names across samples is unsound — a
sample legitimately carries no row for an element it has no insertions of. So each report is
stamped with a digest of its IS table (`# ijump-is-table: <digest>`, a leading comment line
that leaves the table itself byte-for-byte where it was), and `combine-results` compares the
stamps before reading any report as data. New module `report_provenance` owns both ends of
that format.

A report without the stamp is refused rather than merged on trust. That is defensible only
because the two arrive together: any report predating this change also predates the cluster
columns, so its names mean called loci and were never mergeable with names that mean
elements. The message says to rerun the sample.

**Review follow-ups.** Two hard findings, both fixed:

- `combine_results`' lab format still derived its element column by stripping a numeric copy
  suffix — correct for a locus name, ruinous for a cluster name, which has no suffix to
  strip: `ISAba1`, `ISAba11`, `ISAba12` and `ISAba18` all became `ISAba` and their
  frequencies summed into one row. Clustering does that collapsing now, so there is nothing
  left to derive and the column is the reported name. It is `add_element_column`, with a
  regression test, because a comprehension inside `main()` had no seam to test through —
  which is why the change slipped past the first pass.
- The stamp concept had four names across module, constant and prose. It is **annotation
  stamp** now, registered in `CONTEXT.md`, and the module is `annotation_stamp`. The digest
  itself stays `is_table.fingerprint` — the value and the line that carries it are different
  things.

Also: `README` claimed both average-mode files carry the stamp; only the two that
`combine-results` merges do, and `ijump_sum_by_reg.txt` does not. `docs/Precise.md` gained
the same note `docs/Average.md` got. `docs/algorithm.md`'s row describing the copy-suffix
collapse was stale.
