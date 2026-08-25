# 04 — Similarity clusters computed and written to the table

**What to build:** The IS table gains a cluster column grouping loci that are the same
mobile element, computed from sequence similarity rather than from a name prefix.

`isfinder-db-parse` gains a reference-FASTA input, extracts each called locus, runs
all-vs-all BLAST across them, and writes the resulting cluster assignment as a column.
Nothing consumes the column in this ticket — it is verifiable by reading the table.

The rule: single-linkage at **≥95% identity over ≥80% of the shorter locus**, both
exposed as CLI flags. Coverage measured on the *shorter* locus is what lets a short
remnant join its parent element, which is correct — a read clipped at a 76 bp remnant
genuinely cannot be distinguished from one clipped at the full copy. Single linkage is
required, not preferred: on the reference genome the two `IS17` fragments do not align to
each other at all and reach their parent only transitively. Complete linkage would break
the exact case this exists to fix.

Single linkage invites chaining, so every formed cluster is re-checked for internal pairs
that fail the threshold and any such pair is logged by name. Chained merges are visible
rather than silent; the operator edits the cluster column before running the pipeline.

Cluster naming is part of this ticket because a cluster column needs values. A cluster
takes its representative member's base IS name. Suffixes are appended **only on
collision**, ordered by descending cluster size then coordinate. Collisions are real, not
theoretical: the database is sparse relative to actual IS diversity, so two loci ~85%
apart can share a nearest database neighbour and both claim its name. Uniqueness is a
correctness requirement — the coordinate store is a dictionary keyed on this name, so a
duplicate silently overwrites a locus, and the grouping in precise mode would silently
pool two distinct clusters.

**Blocked by:** 03 — IS table carries family and group.

**Status:** ready-for-agent

- [ ] `isfinder-db-parse` accepts a reference FASTA and extracts each called locus from it
- [ ] All-vs-all search over the extracted loci; clusters formed by single linkage at the configured thresholds
- [ ] Identity and coverage thresholds are CLI flags with the stated defaults
- [ ] Each formed cluster is re-checked and every internal pair failing the threshold is logged by locus name
- [ ] Clusters are named for their representative member; suffixes appear only when two clusters would otherwise share a name
- [ ] Cluster names are unique across the table, verified by an explicit test
- [ ] On the reference genome the two `IS17` fragments and `ISAba12_1` share one cluster, and `ISAba53_1` (83% identical) does not join it
- [ ] The three identical `ISAba11` copies share one cluster
- [ ] The new term is added to the domain glossary
