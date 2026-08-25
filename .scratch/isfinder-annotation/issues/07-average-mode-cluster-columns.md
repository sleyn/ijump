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

**Status:** ready-for-agent

- [ ] Per-region reporting emits one column per cluster, named for the cluster
- [ ] The per-region report golden is re-pinned, with the three `IS17`/`ISAba12` columns collapsing to one
- [ ] The multi-sample merge detects reports built from different IS tables and refuses with a clear message
- [ ] Merging reports that do share an IS table behaves as before
- [ ] The output-shape change is documented as breaking
