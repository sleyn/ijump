# Characterization goldens

Committed expectations that pin today's IS-annotation behaviour exactly as it is
(isfinder-annotation 01), so that every later ticket in that series lands as a
reviewable diff under this directory instead of as an assertion someone has to
take on faith.

These are **characterization, not specification**. They record what the code
does, quirks included — not what it ought to do.

## Two tiers

| tier | test | inputs | runs where |
| --- | --- | --- | --- |
| parser-level | `tests/test_golden_isfinder_db_parse.py` | `tests/fixtures/isfinder/` (committed, ~30 KB) | anywhere BLAST+ is installed |
| end-to-end | `tests/test_golden_end_to_end.py` | `Test/Sample.bam` (840 MB) + reference, GFF, IS table | maintainer's machine only |

`tests/fixtures/isfinder/blast.out` is a byte-identical copy of the gitignored
`Test/blast.out` — the ISFinder BLAST outfmt-6 run of `Test/A_baumannii_assembly.fna`
against a clone of the ISFinder database, re-derivable by repeating that search.

`tests/fixtures/isfinder/reference.fna.gz` stands in for that 4 MB assembly, which
clustering needs to extract the called elements from. It carries the same contigs at their
real names and lengths and the called elements at their real coordinates and sequences,
with every other base masked to `N` — so extraction and the all-vs-all search see exactly
what they would on the real assembly, and the parser golden is byte-identical either way.
"Called elements" means **every back-end's**, ISEScan's included: `new_269` has no ISFinder
counterpart, and a span masked to `N` would cluster against nothing.
Rebuild it with `python tests/fixtures/isfinder/make_reference_fixture.py`.

`tests/fixtures/isfinder/legacy_is_table.txt` is the same 13 loci in the legacy four-column
format — a byte copy of the four columns the golden table shares with it. It is what
`tests/test_migrate_is_table.py` migrates, so that the migration back-end is checked against
what the primary back-end produces from the same genome.

`tests/fixtures/isfinder/isfinder_db.fna` is a **stand-in** for the ISFinder database, which
is not ours to redistribute. It carries one entry per element, built from the golden's own
loci and named the way real subject ids are. It stands in for extraction, search, best-hit
selection and subject-id splitting — but not for `pident`, since a locus searched against a
database whose entry *is* that locus matches itself at 100%. Rebuild it with
`python tests/fixtures/isfinder/make_isfinder_db_fixture.py`; the script's docstring has the
full reasoning.

`tests/fixtures/isescan/isescan_results.tsv` is a byte copy of a real ISEScan run over the
same assembly (`Test/ISEScan/`, gitignored). It is what `tests/test_isescan_convert.py`
converts.

The parser tier **skips** where BLAST+ is missing: clustering runs an all-vs-all `blastn`
over the extracted elements. Its inputs are still all committed.

The two tiers are chained: the end-to-end run takes its `--isel` IS table from the
parser tier's own golden, not from a copy in the input directory. Re-pin the parser
golden and the end-to-end tier consumes the new table, instead of silently going on
with a stale one.

The end-to-end tier **skips** — it does not fail, and does not error — wherever
the large inputs or BLAST+ are missing. `tests/test_golden_skip_conditions.py`
covers that skip check itself, and needs none of the large inputs.

Input directory for the end-to-end tier is `Test/` by default; point
`IJUMP_E2E_DATA` at another directory to override. A run takes a couple of
minutes per mode. Deselect the tier with `pytest -m "not e2e"`.

## What is pinned

`e2e/average/` — `ijump_junctions.txt`, `ijump_sum_by_reg.txt`,
`ijump_report_by_is_reg.txt`.

`e2e/precise/` — `ijump_junctions.txt` and `ijump_junction_pairs.txt`. Precise
mode writes no per-region tables.

One file a run produces is deliberately left out: `reads.txt`, a ~900 KB dump
of every clipped read rather than a report table.

## Re-pinning an intended change

```
python tests/regenerate_goldens.py          # both tiers
python tests/regenerate_goldens.py parser   # parser-level only
python tests/regenerate_goldens.py e2e      # end-to-end only
```

Then read `git diff tests/goldens/`. That diff is the deliverable: it is where
an intended change — say the three `IS17_1` / `IS17_2` / `ISAba12_1` columns of
`ijump_report_by_is_reg.txt` collapsing into one — becomes visible instead of
being asserted.
