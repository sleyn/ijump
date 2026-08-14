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
