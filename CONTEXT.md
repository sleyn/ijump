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
`ijump isfinder-db-parse`, read by `ijump run` (`--isel`), and editable by hand in between.
Headered since the annotation columns were added; four-column headerless tables predate that
and are still read. In code it is `is_table.py`, `ISClipped.is_table`, and — for the
coordinate columns alone — `ISClipped.is_coords`.
_Avoid_: mobile elements coordinates file, IS coordinates file, IS list

**Cluster**:
The set of IS table rows that are copies of one mobile element, computed from sequence
similarity rather than from the rows' names — single linkage at ≥95% identity over ≥80% of
the shorter element. A cluster is named for its longest member's base IS name (`ISAba12`),
suffixed `.a`/`.b` only when two clusters would otherwise share a name. It is the `cluster`
column of the IS table, written by `ijump isfinder-db-parse` and editable by hand.
_Avoid_: group (the ISFinder group is separate annotation and never merges anything),
family, IS type, merge key

**Locus**:
One row of the IS table — a single copy of a mobile element as called in the reference,
with its own name, contig and coordinates. Several loci can be one **Cluster**: a full copy
and its own truncated fragments are separate loci of one element.
_Avoid_: hit, IS copy, element (an element may span several loci), region

