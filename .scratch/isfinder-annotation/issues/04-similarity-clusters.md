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

**Status:** done

- [x] `isfinder-db-parse` accepts a reference FASTA and extracts each called locus from it
- [x] All-vs-all search over the extracted loci; clusters formed by single linkage at the configured thresholds
- [x] Identity and coverage thresholds are CLI flags with the stated defaults
- [x] Each formed cluster is re-checked and every internal pair failing the threshold is logged by locus name
- [x] Clusters are named for their representative member; suffixes appear only when two clusters would otherwise share a name
- [x] Cluster names are unique across the table, verified by an explicit test
- [x] On the reference genome the two `IS17` fragments and `ISAba12_1` share one cluster, and `ISAba53_1` (83% identical) does not join it
- [x] The three identical `ISAba11` copies share one cluster
- [x] The new term is added to the domain glossary

## Comments

- Clustering lives in one new module, `src/ijump/is_clustering.py`. Its algorithmic half
  (`cluster_loci`, `name_clusters`) is pure — every alignment it works from is passed in —
  so the linkage and naming rules are tested against the measured all-vs-all numbers the
  spec records, without BLAST or the 4 MB genome. `annotate` is the whole subsystem from a
  back-end's point of view: table plus reference in, `cluster` column out, which is what
  tickets 08 and 09 will reach for.
- `-r/--ref` is required rather than optional. A table with no cluster column is a hard
  error in ticket 06, so a parser that quietly emits one would only be manufacturing that
  error later.
- The parser golden shows the change directly: `IS17_1`, `IS17_2` and `ISAba12_1` all read
  `ISAba12`, `ISAba53_1` stays on its own, the three `ISAba11` copies share `ISAba11`.
- That needed a reference the golden could run against, and `Test/A_baumannii_assembly.fna`
  is gitignored. `tests/fixtures/isfinder/reference.fna.gz` is a stand-in: real contig
  names and lengths, real sequence at the called coordinates, every other base masked to
  `N`, bgzipped to ~20 KB. Verified to produce a byte-identical table to the real assembly.
  `make_reference_fixture.py` beside it rebuilds it.
- The parser tier now needs BLAST+ on PATH and skips without it. Its inputs are still all
  committed; only the all-vs-all search is new.
- Chaining is reported on real data, not just in principle: the `ISAba12` cluster logs
  `IS17_1`/`IS17_2` on every run. That is the *wanted* chain, and nothing in the alignment
  distinguishes it from an unwanted one — hence logged rather than judged.
- The `cluster` column sits between `group` and `pident`, as the spec's column order has it.
  `loci_from_table` rejects a table with two rows of one name, since linkage here and
  `is_coords` later both key on it.
- End-to-end goldens for both modes were re-run against the widened table and are
  byte-identical — nothing consumes the column yet, which is ticket 06 and 07.
