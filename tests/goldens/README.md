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
| parser-level | `tests/test_golden_isfinder_db_parse.py` | `tests/fixtures/isfinder/blast.out` (committed, 5 KB) | anywhere |
| end-to-end | `tests/test_golden_end_to_end.py` | `Test/Sample.bam` (840 MB) + reference, GFF, IS table | maintainer's machine only |

`tests/fixtures/isfinder/blast.out` is a byte-identical copy of the gitignored
`Test/blast.out` — the ISFinder BLAST outfmt-6 run of `Test/A_baumannii_assembly.fna`
against a clone of the ISFinder database, re-derivable by repeating that search.

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
`ijump_report_by_is_reg.txt`, and the five Circos inputs (`karyotype.txt`,
`text.txt`, `links.txt`, `histogram.txt`, `depth.txt`).

`e2e/precise/` — `ijump_junctions.txt` and `ijump_junction_pairs.txt`. Precise
mode writes no per-region tables, and `ijump.py` gates `--circos` on average
mode, so there are no Circos files to pin here.

Two files a run produces are deliberately left out: `reads.txt`, a ~900 KB dump
of every clipped read rather than a report table, and `ijump_data/circos.conf`,
which embeds the absolute path of the run directory and so differs on every run
by construction.

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
