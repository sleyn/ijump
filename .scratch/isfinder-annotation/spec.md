# IS annotation subsystem

Replaces the two ISFinder parser scripts with a single IS-annotation stage: a
defined table contract, pluggable back-ends, and annotation beyond bare
coordinates. Produced by a grilling session on 2026-08-24; every decision below
was made by the maintainer during that session.

## Ticket numbering

**Reference these as `(isfinder-annotation NN)` in commit messages — never bare
`(ticket NN)`.** Numbers `01-09` already exist in `.scratch/isclipped-refactor/`,
`.scratch/packaging/` and `.scratch/review-followups/`.

## Why now

Three problems, each independently confirmed against `Test/` data.

**The HTML path is dead.** `isfinder_parser.py` scrapes an ISFinder BLAST
results page. The site has been down for months, so nobody can produce its
input, and the output had its own defects (see `review-followups 01` and the
`re.DOTALL`-as-`count` bug noted in `isfinder_parser.py` itself).

**Family and group are already in hand and thrown away.** `sseqid` in an
ISFinder BLAST outfmt-6 run is `name_family_group` — `ISAba18_IS3_IS51`,
`ISAba1_IS4_IS10`. `isfinder_db_parcer.py:41` reduces it to the name with
`x.split("_")[0]`. That discards the family and group at zero saving, and
truncates the 11 database entries whose *name* contains an underscore
(`ISBj2_B_IS5_IS5` becomes `ISBj2`).

**The current merge key produces wrong groupings.** Precise mode collapses IS
copies by stripping a `_\d+` suffix (`isclipped.py:487`) — a suffix that
`isfinder_db_parcer.py` invented in the first place. All-vs-all `blastn` on the
13 loci called for `Test/A_baumannii_assembly.fna`:

| pair | identity | aln len | qlen / slen |
| --- | --- | --- | --- |
| `ISAba11_1` ↔ `_2` ↔ `_3` | 100.0% | 1101 | 1101 / 1101 |
| `IS17_2` ↔ `ISAba12_1` | 100.0% | 76 | 76 / 1039 |
| `IS17_1` ↔ `ISAba12_1` | 97.9% | 144 | 144 / 1039 |
| `ISAba53_1` ↔ `ISAba12_1` | 83.0% | 1045 | 1039 / 1039 |
| all other pairs | *no hit* | — | — |

`IS17_1` and `IS17_2` collapse to `IS17` today, but they do not align to each
other at all. In `ISAba12_1`'s frame they are positions 1039→896 and 76→1 — the
two opposite ends of one ISAba12-like element, with the middle 819 bp absent
from the assembly. The correct cluster is `{IS17_1, IS17_2, ISAba12_1}`; the
current key gets both halves wrong.

## Decisions

### Merge key is sequence similarity, not family

Family and group are recorded as annotation and never used to merge.
`ISAba53_1` and `ISAba12_1` share family (IS5) *and* group (IS903) yet are 83%
apart — a 50-150 bp clipped read discriminates them easily, so merging them
would pool two elements the algorithm can currently tell apart. Group is also
unusable as a universal key: 2613 of 5970 database entries have group
`unknown`.

### Clustering

- all-vs-all `blastn` on locus sequences extracted from the reference FASTA.
  BLAST+ is already a hard dependency (`clipped_read_search.py:51`,
  `ijump.py:22`), so this adds nothing to install.
- **≥95% identity over ≥80% of the shorter locus.** Coverage on the *shorter*
  locus is what lets a 76 bp remnant join its parent element — which is correct,
  because a read clipped at that remnant genuinely cannot be told from one
  clipped at the full copy. Both thresholds become CLI flags
  (`--cluster-identity`, `--cluster-coverage`).
- **Single-linkage**, required for fragments to reach a parent they share no
  alignment with (`IS17_1` and `IS17_2` do not align to each other). Complete
  linkage would break the exact case this exists to fix.
- Single-linkage invites chaining: a short fragment landing in a stretch
  conserved between two distinct elements could merge them. So each formed
  cluster is re-checked for internal pairs failing the threshold, and any such
  pair is logged by name. Chained merges are visible, not silent; the operator
  edits the cluster column to override.
- Computed in the parser and written to the table, so it can be inspected and
  hand-edited before `ijump run`.

### Cluster naming

The representative member's base IS name — `{IS17_1, IS17_2, ISAba12_1}` becomes
`ISAba12`. Suffixed `.a`/`.b` (ordered by descending cluster size, then
coordinate) **only on collision**.

Collisions are real, not theoretical. The database is sparse relative to actual
IS diversity, so two loci ~85% apart can share a nearest database neighbour and
both claim its name. `Test/` is one deletion away from showing it: `ISAba53_1`
and `ISAba12_1` stay separate only because ISAba53 happens to be in the
database.

Uniqueness is a correctness requirement, not cosmetics:

- `is_coords` is a dict keyed on IS name (`isclipped.py:296`) — a duplicate key
  silently overwrites the earlier locus.
- Precise mode's `groupby([..., "IS", ...])` would silently pool two clusters,
  which is the conflation this whole change exists to prevent.

### Origin-spanning elements

Detected and flagged (`wraps_origin` plus a shared element id), with both
coordinate rows kept — they are separate spans and the boundary search needs
both. Clustering already merges them by transitivity, so the flag is for
visibility: it tells the reader the assembler broke the contig inside an IS
copy.

Deliberately *not* merged into one `start > stop` row. `set_is_boundaries`,
`circos.py`'s span drawing and `region_summary`'s overlap logic all assume
`start <= stop`, and a joined row would claim a 220 bp element where the truth
is 1039 bp with a hole. Note `isclipped.py:328` already assumes contigs are
circular (`ASSUMPTION OF COMPLETENESS!`) for radius windows, though nothing in a
SPAdes FASTA states it.

### Table format

Headered TSV: `is_name, contig, start, stop, family, group, cluster, pident,
wraps_origin, element_id`. Legacy headerless 4-column files are sniffed and
accepted by `iscollect`.

`circos.py:102` does `is_chrom, is_start, is_stop = is_coords[is_name]` — a
strict 3-tuple unpack that any added column breaks. Fixing it is part of the
format change, not a follow-up.

### Both modes group on cluster

Precise mode replaces the `_\d+` strip at `isclipped.py:487`. Average mode emits
one `region_summary` column per cluster instead of one per IS name — today it
gives `Test/` separate `IS17_1`, `IS17_2` and `ISAba12_1` columns for one
element and two of its own fragments.

**This changes the shape of `ijump_report_by_is_reg.txt`**, a documented
collectable output (`CONTEXT.md`). Breaking for anyone with downstream scripts.

A table with no cluster column is a hard error naming the migration subcommand.
That is acceptable only because a remedy exists; without one it would be
user-hostile.

### Back-ends

- **ISFinder BLAST** (existing, extended): outfmt-6 + reference FASTA → table.
- **Migration**: old table + reference FASTA + ISFinder BLAST database → new
  table. Re-annotates family from a fresh locus-vs-database BLAST, because a
  4-column file carries no family to recover. Preserves hand-curated
  coordinates.
- **ISEScan converter**: reads ISEScan's `.tsv` (byte-identical in content to
  its `.csv`, and the machine-readable form of `.out`). iJump never invokes
  ISEScan.

All three share one annotate-and-cluster core.

### ISEScan is read, never run

The two annotations are complementary on `Test/A_baumannii_assembly.fna`:

- ISEScan finds 3 × `new_269` elements (2158, 2315, 5404 bp) with **no ISFinder
  BLAST hit at all**. An IS absent from ISFinder that is actively jumping is
  invisible to iJump today. That is the sensitivity argument.
- ISFinder BLAST finds both `IS17` contig-edge fragments. ISEScan needs TIRs
  plus an ORF, so a 76 bp remnant is structurally invisible to it.
- Where both fire they disagree on span: `ISVsa3_1` 110032-111008 (977 bp) vs
  `IS91_50` 109036-111334 (2299 bp), same locus. Span drives
  `set_is_boundaries`' search windows, so this is not cosmetic.

Depending on ISEScan was rejected on cost: x64-only (needs Rosetta on arm64), a
`libgsl.25` symlink workaround, ~14 min on this genome — in a project whose
`CLAUDE.md` already documents that `uv sync` cannot resolve the runtime deps at
all. Reading its output puts that burden only on operators who choose it.

A **union mode** merging both back-ends would be strictly the best annotation
for this genome, but needs a documented rule for the ISVsa3/IS91_50 span
conflict — a scientific judgement, not a coding one. Deferred, not rejected.

## Verification

Characterization goldens are pinned **before** any change, so each commit shows
a reviewable diff and the intended changes (three `IS17`/`ISAba12` columns
collapsing to one) appear as a diff rather than an assertion.

Baseline confirmed reproducible on 2026-08-24: re-running `blastn` against a
fresh `Test/ISfinder-sequences` clone and re-running `isfinder_db_parcer.py`
reproduces `Test/ISTable_processing.txt` byte-for-byte. The BLAST output differs
from `Test/blast.out` only in `pident` formatting (`100` vs `100.000`); same 49
records.

`Test/Sample.bam` is 840 MB, so end-to-end goldens cannot live in CI and are
local-only. Parser-level goldens (`blast.out` → IS table) are small and
committable. CI does not run pytest at all today (`CLAUDE.md`).

## Out of scope

**The 75% overlap rule in `isfinder_db_parcer.py` is unchanged.** It is what
lets `ISAlw32_1` (533 bp of a ≥916 bp element) and `ISAcsp3_1` (734 bp of a ≥720
bp one) sit adjacent as two truncated calls. Truncated boundaries feed
`set_is_boundaries` directly, so this is a real accuracy question — worth its
own investigation, not a rider on this work.

**Union of ISFinder and ISEScan back-ends.** See above.

## Open risk

`combine_results.py` merges per-sample reports on IS identity. That is safe only
if every sample in a combine shares one reference and one IS table. Combining
samples annotated against different references would leave cluster names
unaligned — worse than today, because cluster names are derived rather than
fixed. Needs a check where average mode switches to cluster columns.
