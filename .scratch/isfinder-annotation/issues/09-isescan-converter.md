# 09 — Convert ISEScan output to an IS table

**What to build:** A subcommand reading ISEScan's tab-separated results, plus a reference
FASTA, into the same IS-table contract — annotated and clustered by the same rules. iJump
never invokes ISEScan; the operator runs it and hands over the output.

The two annotations are complementary on the reference genome, which is the whole argument
for supporting it: ISEScan finds three copies of an element with **no ISFinder database hit
at all**, so an IS that is absent from ISFinder and actively jumping is invisible to iJump
today. In the other direction ISEScan needs terminal repeats plus an ORF, so the short
`IS17` remnant is structurally invisible to it.

Reading rather than running was a cost decision: ISEScan is x64-only (needs emulation on
arm64), wants a library symlink workaround, and takes ~14 minutes on this genome — in a
project whose runtime dependencies already cannot be resolved by the standard installer.
Reading its output puts that burden only on operators who choose it.

Elements ISEScan names but ISFinder does not know keep their ISEScan family rather than
being dropped.

Note for the record: where both back-ends fire they disagree on span — one locus is called
977 bp by one and 2299 bp by the other. Span drives the boundary search windows, so this is
not cosmetic. Unioning the two back-ends is explicitly **out of scope** here; it needs a
documented rule for that conflict, which is a scientific judgement.

**Blocked by:** 04 — Similarity clusters computed and written to the table.

**Status:** ready-for-agent

- [ ] A subcommand reads ISEScan's tab-separated output and a reference FASTA, emitting a new-format table
- [ ] Elements with no ISFinder counterpart are retained, carrying their ISEScan family
- [ ] Clusters are assigned by the same shared core and thresholds as the other back-ends
- [ ] The resulting table runs through both estimation modes unmodified
- [ ] Documentation states that iJump reads ISEScan output and never invokes ISEScan
