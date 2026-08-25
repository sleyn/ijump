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

**Status:** done

- [x] The subcommand is named `migrate-is-table` — 06's error already points at it
- [x] A subcommand accepts an old-format table, a reference FASTA and the ISFinder database, and emits a new-format table
- [x] Input coordinates and contig assignments are carried through byte-identically
- [x] Family and group are re-derived from a fresh locus-versus-database search
- [x] Clusters are assigned by the same rules and thresholds as the primary back-end
- [x] Migrating the committed legacy table reproduces what the primary back-end produces from the same genome
- [x] The annotate-and-cluster logic is shared with the primary back-end, not duplicated

## Comments

**This unblocks the release.** Tickets 06 and 07 make `ijump run` refuse a table with no
cluster column, naming `ijump migrate-is-table` as the remedy. The remedy now exists, so
that error is no longer a dead end.

**The shared core is `is_annotation.annotate_and_cluster`**, called by both back-ends. Two
smaller things that were about to be duplicated went with it: the argparse threshold type
(`is_clustering.threshold_type`, which rejects `--cluster-coverage 80` rather than silently
answering wrongly) and the BLAST runner, promoted from `is_clustering._run` to
`run_blast_command` now that a second module searches.

**Reproduction is verified on every column but `pident`.** The committed legacy table is the
same 13 loci as the golden, and migrating it reproduces the golden's `family`, `group`,
`cluster`, `wraps_origin` and `element_id` exactly. `pident` is excluded on purpose: the
tests search a **stand-in** ISFinder database, since the real one is not ours to
redistribute, and a locus searched against a database whose entry *is* that locus matches
itself at 100% where the real database gives 98.556. `tests/fixtures/isfinder/README` — the
generator's docstring — states what the stand-in can and cannot stand for. Reproducing
`pident` needs the real database.

Worth noting for whoever gets the real database: byte-identical `pident` may not be
achievable even then. The golden's value comes from searching the *genome* against ISFinder,
where the alignment has genomic context either side of the locus; migration searches the
extracted locus alone. The two are the same search of different queries.

**A locus that matches nothing is kept**, with family, group and pident empty and a warning
naming it. This back-end annotates; it does not re-call loci, and dropping a row the
operator curated by hand would be exactly the re-calling it promises not to do. Clustering
still groups it, since clustering reads the sequences rather than the database.

**Fixture note.** This build of `makeblastdb` rejects a small synthetic database (a single
random sequence, wrapped or not, in any directory), while accepting the 10-entry fixture
built from real sequence. The no-hit test therefore stubs the search rather than pointing at
a database guaranteed to miss — what is under test is what `migrate` does with a miss, and
`_best_hits` is unit-tested separately against a BLAST output file with a miss in it.
