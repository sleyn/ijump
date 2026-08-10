# 02 — Stop terminating the process from inside the library

Status: ready-for-agent
Blocked by: 01

## Why

`isclipped.py` calls `exit()` at four places. Three of them guard the no-results paths, so
no test can reach past them — the process dies mid-assertion. Nine of the last twenty-five
commits fix empty-input handling; the same class of bug keeps returning because those paths
are unverifiable.

It also carries a user-visible bug. `check_data_presence_in_df` (ijump.py:28) is called in
the common prefix at ijump.py:189, **before** the mode branch at ijump.py:201, and always
calls `empty_pairs_out()` (ijump.py:22). So an average-mode run with no clipped reads emits
the precise-mode pairs schema and never writes `ijump_report_by_is_reg.txt` — the file
`README.md:244` instructs users to collect for `combine_results.py`. Such a sample silently
drops out of the combined table.

Architecture review candidate 01. Vocabulary: `/codebase-design`.

## Scope

### 1. Introduce the exception

Define in `isclipped.py`, next to `ISClipped`:

```python
class NoInsertionsFound(Exception):
    """The analysis ran correctly and there is nothing to report."""
```

No separate exceptions module for a single exception. `ijump.py` already imports from
`isclipped`.

### 2. Convert all six termination sites

| Site | Today | After |
|---|---|---|
| isclipped.py:490 | `exit(0)` no BLAST output file | `raise NoInsertionsFound` |
| isclipped.py:539 | `exit(0)` no significant BLAST hits | `raise NoInsertionsFound` |
| isclipped.py:581 | `exit(0)` no hits outside IS boundaries | `raise NoInsertionsFound` |
| ijump.py:32 | `exit(0)` no clipped reads | `raise NoInsertionsFound` |
| ijump.py:46 | `exit(0)` no junctions | `raise NoInsertionsFound` |
| isclipped.py:892 | `exit(1)` bad `orientation` argument | `raise ValueError` |

`isclipped.py:892` is a **programmer-error guard, not an empty-data outcome**.
`_read_count_mtx` is only ever called with the literals `'left'` and `'right'`
(isclipped.py:1036, :1039), so it is unreachable. It must NOT raise `NoInsertionsFound` —
reporting a caller bug as "no insertions found" with exit 0 would turn a crash into silent
wrong output.

Remove the inline `pairs_table_empty().to_csv(...)` writes at isclipped.py:488-489,
:537-538 and :579-580. The library no longer writes output files.

### 3. One catch site in the driver

```python
def main():
    ...
    try:
        run_pipeline(args)
    except NoInsertionsFound as e:
        logging.info(str(e))
        write_empty_outputs(args.estimation_mode, output_dir, is_processing)
        sys.exit(0)
```

Exit code stays **0** — a run that finds nothing is a successful run (see `CONTEXT.md`).

### 4. The file-set invariant

> A run that finds nothing writes the same files as a run that finds something, with
> headers and zero data rows.

| Mode | Files written |
|---|---|
| average | `reads.txt`, `ijump_junctions.txt`, `ijump_sum_by_reg.txt`, `ijump_report_by_is_reg.txt` |
| precise | the above, plus `ijump_junction_pairs.txt` |

All five schemas already exist — no new schema code:

| File | Initializer |
|---|---|
| `reads.txt` | `ISClipped._cltbl_init()` (isclipped.py:182) |
| `ijump_junctions.txt` | `ISClipped._jtbl_init()` (isclipped.py:228) |
| `ijump_sum_by_reg.txt` | `ISClipped.sum_by_reg_tbl_init()` (isclipped.py:213) |
| `ijump_report_by_is_reg.txt` | `ISClipped.report_table_init()` (isclipped.py:171) |
| `ijump_junction_pairs.txt` | `ISClipped.pairs_table_empty()` (isclipped.py:137) |

`sum_by_reg_tbl_init` and `report_table_init` are currently called at ijump.py:209 and
ijump.py:216 with their return values discarded. This gives them their first real use;
clean up the discarded calls.

`write_empty_outputs` lives in `ijump.py` — the driver owns I/O. Do not create an outputs
module; that is architecture review candidate 03.

## Out of scope

- Any change to `ISClipped`'s interface or call order (candidate 02).
- A shared schema module (candidate 03).
- The `'presice'` typo at ijump.py:39 — ticket 03. **Do not fix it here.** It would shift
  `IS pos` values in `ijump_junctions.txt` for precise-mode runs, contaminating the
  before/after baseline that makes this refactor verifiable.

## Verification

Uses `tests/` and `FakeAlignment` from ticket 01.

**Seam tests** — `tests/test_no_results_paths.py`, one per site, using `FakeAlignment`.
No BAM, no GFF, no reference:

```python
def test_missing_blast_output_signals_no_insertions(tmp_path):
    isc = ISClipped(FakeAlignment(), 'ref', 'unused.gff', str(tmp_path))
    with pytest.raises(NoInsertionsFound):
        isc.parseblast('missing.out', 1)
```

Cover isclipped.py:490 (missing/empty BLAST output), :539 (all hits below
`blast_min_ident`, isclipped.py:506), :581 (`make_gene_side_regions` with all hits inside
IS boundaries), and `_read_count_mtx` raising `ValueError` on a bad `orientation`.

**End-to-end** — invert the three "does not exist" assertions ticket 01 pinned, and assert
the full invariant: same file set in both modes, every file present, every file 0 data rows,
exit code 0.

**Baseline** — on real data with results, output must be byte-identical before and after.
Only the empty-run file set changes.

## Done when

- No `exit()` remains in `isclipped.py`.
- Both modes produce their full file set on an empty run.
- `pytest` passes with no BLAST+ installed.

## Comments

Implemented in bf743ef (`fix: Stop terminating the process from inside the library`): removed
the `exit()` call from `isclipped.py` and wired `write_empty_outputs` so both estimation modes
produce their full file set on an empty run.
