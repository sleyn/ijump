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

**Status:** done

- [x] A subcommand reads ISEScan's tab-separated output and a reference FASTA, emitting a new-format table
- [x] Elements with no ISFinder counterpart are retained, carrying their ISEScan family
- [x] Clusters are assigned by the same shared core and thresholds as the other back-ends
- [x] The resulting table runs through both estimation modes unmodified
- [x] Documentation states that iJump reads ISEScan output and never invokes ISEScan

## Comments

**The subcommand is `ijump isescan-convert`.** It reads ISEScan's `.tsv` and a reference
FASTA; everything after the four locus columns is `is_annotation.annotate_and_cluster`, the
same core the other two back-ends call. That is now checked by watching the call rather than
by re-deriving the output, so the test fails if this back-end ever grows its own copy of the
rules.

**ISEScan's `cluster` column and iJump's are different things, and the converter takes
iJump's.** This surfaced as a failing test I had written on the assumption they agree. On
the reference genome ISEScan files all three `new_269` calls under one of its clusters;
read off the sequences they are two by iJump's rule:

- `new_269_1` (2158 bp) aligns to `new_269_3` (5404 bp) over its whole length at **90.7%**
  identity — under the 95% that says a clipped read could not tell them apart.
- `new_269_1` aligns to `new_269_2` (2315 bp) not at all.
- `new_269_2` and `new_269_3` do meet the rule, at 96.4% over 88% of the shorter.

So the table carries `new_269.a` and `new_269.b`. That is the documented rule working, not a
disagreement to paper over, and an operator who trusts ISEScan can edit the column. Both the
test and the README say so.

**Naming.** ISEScan reports no element name, only a family and its own cluster id, so a
locus is named `<ISEScan cluster>_<copy number>` — the `<element>_<n>` shape the other
back-ends use, so `base_is_name` strips the counter and leaves the element. `group` and
`pident` are left empty: ISEScan reports neither an ISFinder group nor an identity against
any database, and filing its numeric cluster id under `group` would put a different kind of
value under a name that already means something.

**"Runs through both estimation modes" is verified by doing it.** `tests/test_isescan_end_to_end.py`
converts the committed ISEScan results against the real assembly and drives both modes on
the sample alignment. Nothing is pinned: the ISEScan annotation calls different spans by
design, so its outputs are not the goldens, and pinning them would be pinning a second
baseline nobody reads. `run_e2e_pipeline` grew an optional IS-table argument for it; the
goldens still pass nothing and still chain to the parser tier.

**Two fixture repairs came with this.**

- `reference.fna.gz` masks everything outside a called locus to `N`, and "called" meant the
  ISFinder parser's loci alone. `new_269` has no ISFinder counterpart, so its spans were all
  `N` and would have clustered against nothing. The generator now unmasks every back-end's
  loci. The parser golden is byte-identical either way, since the ISFinder spans are
  untouched.
- `make_reference_fixture.py` computed `REPO_ROOT` with two `dirname`s where the file's depth
  needs three, so every path it built pointed inside `tests/` and the script could not run at
  all. It lost a level when the fixtures moved into `isfinder/`. Fixed; the rebuilt fixture
  is in this commit.

**Union is still out of scope**, as filed. Where both back-ends fire they disagree on span —
977 bp against 2299 bp for one locus — and span drives the boundary search windows. The rule
for resolving that is a scientific judgement.

**Review follow-ups.** One hard finding and two duplications that a third back-end tipped
over:

- `README`'s IS-table section still said the annotation columns "are the ISFinder family and
  group" and that "all of them are filled in by the **isfinder-db-parse** subcommand".
  Both halves are wrong for a table `isescan-convert` wrote. It now names all three
  back-ends, says which columns they agree on (`cluster`, `wraps_origin`, `element_id`) and
  which they do not, and that only `cluster` is required.
- Three back-ends each carried a byte-identical eleven-line `--cluster-identity` /
  `--cluster-coverage` block, so changing a default was an edit in three files. They are
  arguments to the shared core rather than to any back-end, so they are added by
  `is_annotation.add_cluster_arguments` from one place.
- `read_isescan` hand-listed five empty columns to satisfy the reindex, three of which
  `annotate_and_cluster` immediately overwrites — a second copy of a fill rule
  `read_is_table` already had. It is `is_table.with_all_columns` now, used by both, so a
  back-end states only the columns it can speak to.

**Spec-review follow-ups.** Two test gaps closed:

- The both-modes test asserted only exit 0 and three filenames, which a header-only report
  would also pass — leaving the ticket's own sensitivity argument untested. It now asserts
  that `new_269`, the element only this back-end can see, actually reaches the junctions
  report. It does, in both modes.
- Nothing ran `ijump isescan-convert` through the CLI and read the table back off disk, as
  the other back-ends' tests do. Added.

The review also noted, from the fixture, that ISEScan itself reports `ncopy4is=2` for
`new_269` with two of the three rows typed `p` (partial) — so ISEScan does not claim three
complete copies either, and the two-cluster split agrees with its own output better than the
ticket's "three copies" phrasing suggested.
