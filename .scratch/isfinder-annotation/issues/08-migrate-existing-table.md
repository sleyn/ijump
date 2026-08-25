# 08 — Migrate an existing IS table to the new format

**What to build:** A subcommand that takes an IS table in the old format, a reference FASTA
and the ISFinder BLAST database, and produces a table in the new format — so an operator
with a hand-curated table is not forced to regenerate it from scratch.

Hand-curated coordinates are preserved exactly; the subcommand adds annotation, it does not
re-call loci. Family and group are re-derived by searching each locus against the ISFinder
database, because a 4-column file carries no family to recover. Clustering runs over the
loci by the same rules as the primary back-end.

This is the remedy that the hard error in precise mode points at, so it must be usable
before that error can ship.

All three back-ends — this one, the ISFinder BLAST parser, and the ISEScan converter —
share one annotate-and-cluster core rather than reimplementing the rules three times.

**Blocked by:** 04 — Similarity clusters computed and written to the table.

**Status:** ready-for-agent

- [ ] A subcommand accepts an old-format table, a reference FASTA and the ISFinder database, and emits a new-format table
- [ ] Input coordinates and contig assignments are carried through byte-identically
- [ ] Family and group are re-derived from a fresh locus-versus-database search
- [ ] Clusters are assigned by the same rules and thresholds as the primary back-end
- [ ] Migrating the committed legacy table reproduces what the primary back-end produces from the same genome
- [ ] The annotate-and-cluster logic is shared with the primary back-end, not duplicated
