# iJump

iJump finds insertion-sequence (IS) rearrangements in evolved bacterial populations by
examining clipped reads in a short-read alignment. This glossary is seeded lazily — terms
appear here as they get resolved, not all at once.

## Language

**No insertions found**:
The outcome where the analysis ran correctly and there is nothing to report — a sample
genuinely carrying no detectable IS rearrangements. Distinct from a failure. A run that
finds nothing is a successful run: it exits 0 and writes the same output files as a run
that finds something, with zero data rows.
_Avoid_: empty result, no data, failure, error

**Average mode**:
The pipeline variant that estimates IS frequency by averaging region coverage against the
number of clipped reads, reporting one row per IS-and-genetic-element pair. Its collectable
output is `ijump_report_by_is_reg.txt`.
_Avoid_: approximate mode, fast mode

**Precise mode**:
The pipeline variant that attempts to resolve each insertion event separately by pairing
left and right junctions. Its collectable output is `ijump_junction_pairs.txt`.
_Avoid_: exact mode, detailed mode

**Region**:
An annotated genetic element — a gene or an intergenic interval — that junctions get grouped
into for average-mode reporting. One row of `ijump_sum_by_reg.txt`/`ijump_report_by_is_reg.txt`
is one region.
_Avoid_: gene (a region may be intergenic), locus, GE

**IS table**:
The tab-separated file listing the IS elements a run works from — one row per copy found in
the reference, carrying its name, contig, coordinates and ISFinder annotation. Written by
any IS-table back-end (`ijump isfinder-db-parse`, `ijump migrate-is-table`,
`ijump isescan-convert`), read by
`ijump run` (`--isel`), and editable by hand in between.
Headered since the annotation columns were added; four-column headerless tables predate that
and are still read. In code it is `is_table.py`, `ISClipped.is_table`, and — for the
coordinate columns alone — `ISClipped.is_coords`.
_Avoid_: mobile elements coordinates file, IS coordinates file, IS list

**Cluster**:
The set of IS table rows that are copies of one mobile element, computed from sequence
similarity rather than from the rows' names — single linkage at ≥95% identity over ≥80% of
the shorter element. A cluster is named for its longest member's base IS name (`ISAba12`),
suffixed `.a`/`.b` only when two clusters would otherwise share a name. It is the `cluster`
column of the IS table, written by every back-end that produces one and editable by hand.
_Avoid_: group (the ISFinder group is separate annotation and never merges anything),
family, IS type, merge key

**Locus**:
One row of the IS table — a single copy of a mobile element as called in the reference,
with its own name, contig and coordinates. Several loci can be one **Cluster**: a full copy
and its own truncated fragments are separate loci of one element.
_Avoid_: hit, IS copy, element (an element may span several loci), region

**Annotation stamp**:
The line at the head of a report naming the IS table the run was annotated against
(`# ijump-is-table: <digest>`), written by `ijump run` and read by `ijump combine-results`.
Cluster names are derived from the loci rather than fixed, so two runs annotated against
different tables can use one name for different elements; the stamp is what lets a merge
tell that its samples share one **Cluster** vocabulary before joining them on it. The digest
is `is_table.fingerprint`; the format lives in `annotation_stamp`.
_Avoid_: provenance, checksum, watermark, version

**Origin-spanning element**:
One copy of a mobile element called as two **Locus** rows because the assembler broke its
contig in the middle of it — one row ending at the contig's last base, one starting at its
first. Both rows are kept as they are; they carry `wraps_origin` and a shared `element_id`
in the **IS table**, which say the boundary is the assembly's, not the genome's. Detection
needs the two rows to share a **Cluster** and a contig and to sit at its opposite ends.
_Avoid_: wrapped element, circular junction, split IS, broken contig
