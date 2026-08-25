# 05 — Origin-spanning elements detected and flagged

**What to build:** When one IS copy is split across a contig boundary — one span ending at
a contig end, another starting at that contig's origin, both in the same cluster — the
table says so: a flag on each row plus a shared element id linking the two.

Both coordinate rows are kept. They are genuinely separate spans and the boundary search
needs both. Deliberately *not* merged into a single row with start greater than stop: the
boundary search, the Circos span drawing and the per-region overlap logic all assume start
precedes stop, and a joined row would claim a 220 bp element where the truth is a 1039 bp
element with a hole in the middle.

Clustering already unites these loci by transitivity, so the flag adds no grouping — its
job is visibility. It tells whoever reads the table that the assembler broke a contig
inside an IS copy, which is information currently invisible everywhere in the output.

**Blocked by:** 04 — Similarity clusters computed and written to the table.

**Status:** ready-for-agent

- [ ] Loci in a shared cluster that abut opposite ends of the same contig are detected
- [ ] Each such row carries an origin-spanning flag and an element id shared with its counterpart
- [ ] Both coordinate rows are preserved unchanged; no row is emitted with start greater than stop
- [ ] Loci that merely sit near a contig end without a counterpart at the origin are not flagged
- [ ] On the reference genome the two `IS17` fragments are flagged and share one element id
- [ ] The new term is added to the domain glossary
