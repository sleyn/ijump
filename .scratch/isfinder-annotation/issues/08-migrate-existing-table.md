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

**Reproduction is verified on five of the six annotation columns.** The committed legacy
table is the same 13 loci as the golden, and migrating it reproduces the golden's `family`,
`group`, `cluster`, `wraps_origin` and `element_id` exactly. `pident` is not compared,
because the tests search a **stand-in** ISFinder database — the real one is not ours to
redistribute — and a locus searched against a database whose entry *is* that locus matches
itself at 100% where the real database gives 98.556. `tests/goldens/README.md` says what the
stand-in is; `make_isfinder_db_fixture.py`'s docstring what it can and cannot stand for.

**Be precise about what the family/group check proves.** The stand-in's entries are built
from these loci and labelled from this golden, so the check proves the plumbing — extract,
search, best hit, split the subject id — and that each locus picks an entry carrying its own
family and group. It cannot show the labels are biologically right, and it cannot tell
`IS17`, `ISAba12` and `ISAba53` apart, since all three carry `IS5`/`IS903`. Fidelity needs
the real database.

**A `pident` divergence that had nothing to do with the stand-in**, found in review and
fixed: the primary back-end reads its hits through `pd.read_csv`, so its identity reaches
the table as a float and is written `100.0`, while migration kept BLAST's raw `100.000`.
Same number, two spellings in one file depending on which back-end wrote it, on any
database. It goes through `float` now.

**An earlier note here claimed byte-identical `pident` may be unachievable even with the
real database, because the golden's search carries genomic context the extracted locus does
not. That was wrong** — the golden's `start`/`stop` *are* that search's HSP boundaries, so
the extracted locus is the alignment, and blastn is local. The real difference is
hit-selection: the parser filters at `evalue <= 1e-30` with an overlap-dedup pass, migration
takes the best bitscore at 1e-5. That is deliberate and now documented on `search_database`
— the parser is *calling* loci from a whole genome, where a lax threshold invents them,
while migration annotates loci already decided on, the shortest of which is 76 bp and cannot
reach 1e-30 however good the match.


**A locus that matches nothing is kept**, with family, group and pident empty and a warning
naming it. This back-end annotates; it does not re-call loci, and dropping a row the
operator curated by hand would be exactly the re-calling it promises not to do. Clustering
still groups it, since clustering reads the sequences rather than the database.

**Fixture note.** This build of `makeblastdb` rejects a small synthetic database (a single
random sequence, wrapped or not, in any directory), while accepting the 10-entry fixture
built from real sequence. The no-hit test therefore stubs the search rather than pointing at
a database guaranteed to miss — what is under test is what `migrate` does with a miss, and
`_best_hits` is unit-tested separately against a BLAST output file with a miss in it.
