# 01 — Test scaffolding and a characterization test for the empty run

Status: ready-for-agent
Blocked by: —

## Why

`ijump` has no runnable test suite. The only two test files sit in the gitignored 8.7 GB
`Test/` directory: `Test/test_find_pair.py` asserts nothing (it prints `OK` unconditionally
at line 174) and `Test/iJump_test/test_pysam_pileups.py` is 36 lines of scratch. Ticket 02
changes control flow across six sites with no safety net. This ticket builds the net.

See `.scratch/isclipped-refactor/notes/no-test-safety-net.md`.

## Scope

Create a tracked `tests/` directory at the repo root, separate from the gitignored `Test/`.

```
tests/
  conftest.py
  fake_alignment.py
  test_empty_run_outputs.py
  fixtures/
    tiny.fna          few KB, single short contig
    tiny.nsq          BLAST db for tiny.fna, generated once with makeblastdb
    tiny.nhr
    tiny.nin
    tiny.bam          reads that align but produce NO clipped reads
    tiny.bam.bai
    tiny.gff          annotation for tiny.fna
    is_coords.txt     one IS element, tab-separated: name chrom start stop
```

Add pytest configuration (`pyproject.toml` or `pytest.ini`) with `testpaths = tests`.

### FakeAlignment

`ISClipped.__init__` (isclipped.py:21-119) touches only `aln.references` and `aln.lengths`
(isclipped.py:80-83), and `gff.gff(...)` does not read from disk until `readgff()`
(gff.py:5-9, gff.py:12). So constructing `ISClipped` in a test needs no BAM, no reference
and no GFF:

```python
class FakeAlignment:
    references = ('contig_1',)
    lengths = (10000,)
```

Ticket 02 uses this for its seam tests. Provide it here.

### The characterization test

`tests/test_empty_run_outputs.py` runs the full CLI on the fixture in both modes and pins
**what iJump does today**, including the bug. Do not fix anything here.

Today, with no clipped reads, execution reaches `check_data_presence_in_df`
(ijump.py:28, called at ijump.py:189) and exits 0 after writing only
`ijump_junction_pairs.txt` via `empty_pairs_out()` (ijump.py:22) — in **both** modes.

So the test asserts, for both `--estimation_mode average` and `--estimation_mode precise`:

- exit code is 0
- `ijump_junction_pairs.txt` exists and has 0 data rows
- `ijump_report_by_is_reg.txt` does **not** exist
- `ijump_sum_by_reg.txt` does **not** exist
- `ijump_junctions.txt` does **not** exist

Mark the three "does not exist" assertions clearly as pinning current, incorrect behaviour —
ticket 02 inverts them.

## Constraints

- The fixture must produce **no clipped reads**, so execution exits at ijump.py:189 before
  `runblast` (ijump.py:194). With `tiny.nsq` committed, `check_blast_ref` (ijump.py:50)
  short-circuits at line 51 and never spawns `makeblastdb`. Net effect: **the suite runs on
  a fresh clone with no BLAST+ installed.** Verify this by running the suite with BLAST+
  off `PATH`.
- Do not add `.scratch/` or `tests/` to `.gitignore`. `Test/` stays gitignored.
- Total committed fixture size should stay under ~50 KB.

## Done when

- `pytest` passes from a clean clone with no bioinformatics tooling installed.
- The characterization test documents today's behaviour in both modes.
- `FakeAlignment` is importable by ticket 02.

## Comments

Implemented in dabecf6 (`test: Add tracked test scaffolding and empty-run characterization
test`): committed the tiny fixture set, added `FakeAlignment`, and pinned today's empty-run
behaviour in both modes.