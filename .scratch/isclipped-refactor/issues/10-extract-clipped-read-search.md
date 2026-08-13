# 10 — Extract clipped-read search into a standalone `clipped_read_search` module

Status: ready-for-agent
Blocked by: —

## Why

The clipped-read-search cluster — `_clboundaries` (`isclipped.py:325-368`),
`collect_clipped_reads` (`isclipped.py:409-428`), `_cl_read_cov_overlap`
(`isclipped.py:431-447`), `_crtable_ungapped` (`isclipped.py:454-522`),
`_write_cl_fasta` (`isclipped.py:526-538`), `runblast` (`isclipped.py:542-551`),
`_choosecoord` (`isclipped.py:554-560`), `parseblast` (`isclipped.py:564-621`) — is the
detection pipeline's first algorithmic stage: turn an alignment into a BLAST hit table of
candidate junction positions. It reaches into eight pieces of `ISClipped` instance state
(`clipped_reads_dict`, `clipped_reads_bwrd_dict`, `cl_read_cov_overlap`, `match_lengths`,
`read_lengths`, `n_reads_analyzed`, `_index`, `blastout_filtered`) and is threaded together by
a bare `direction` int (1 = IS→Ref forward search, 0 = Ref→IS backward search, precise mode
only). Nothing can exercise any one of these methods without constructing a full `ISClipped`
backed by a real or faked `pysam.AlignmentFile` — the same shape of gap tickets 06/08/09 closed
for `_find_pair` and `assess_isel_freq`.

Architecture review candidate 01 (`/improve-codebase-architecture`, 2026-08-11), grilled to a
shared understanding the same session. Vocabulary: `/codebase-design`.

## Scope

### 1. Create the module

Create `clipped_read_search.py` at the repo root (alongside `junction_pairing.py`,
`frequency_estimation.py`, `gff.py`) with:

```python
@dataclass
class SearchResult:
    clipped_reads: pd.DataFrame       # replaces self.clipped_reads / self.clipped_reads_bwrd
    blast_hits: pd.DataFrame          # replaces self.blastout_filtered
    match_lengths: list[int]          # forward (direction=1) only; [] on backward
    read_lengths: int                 # forward only; 0 on backward
    n_reads_analyzed: int             # forward only; 0 on backward
    cl_read_cov_overlap: dict         # backward (direction=0) only; {} on forward

def search(aln, boundaries, ref_name, workdir, direction, run_blast=run_blast_subprocess) -> SearchResult:
    ...
```

Move the eight methods listed above into this module as module-private functions (drop
`self`, `_clboundaries`/`_choosecoord` stay `@staticmethod`-shaped free functions). Preserve
each body verbatim — relocation, not a rewrite. **Keep the `direction` flag as-is** — do not
split forward/backward into separate functions in this ticket (that split, if wanted, is a
separate future ticket; this one is a mechanical move only).

`run_blast_subprocess` is the one injectable seam: a small function wrapping today's
`NcbiblastnCommandline` invocation (`isclipped.py:549-551`), taking the same file paths
`runblast` does today and returning nothing (it writes `out_file` as a side effect, same as
today). `search`'s `run_blast` parameter defaults to it, so production callers change nothing,
but tests can substitute a fake that writes a canned BLAST output file without invoking
`blastn`.

Local dict-building state (`clipped_reads_dict`/`clipped_reads_bwrd_dict`, `_index`) becomes
local variables inside `search`, not fields threaded from outside.

### 2. Update the call site

`ISClipped.run()` (`isclipped.py:834-939`) calls `clipped_read_search.search(...)` once per
direction (mirrors today's `collect_clipped_reads()` / `runblast()` / `parseblast()` call
sequence at `isclipped.py:839-847` for forward, `isclipped.py:875-883` for backward) and
assigns the returned `SearchResult` fields onto `self.*` at the call site:

```python
result = clipped_read_search.search(self.aln, self.boundaries, self.ref_name, self.workdir, direction=1)
self.clipped_reads = result.clipped_reads
self.blastout_filtered = result.blast_hits
self.match_lengths = result.match_lengths
self.read_lengths = result.read_lengths
self.n_reads_analyzed = result.n_reads_analyzed
```

(backward call sets `self.clipped_reads_bwrd`, `self.cl_read_cov_overlap` instead). Add
`import clipped_read_search` to `isclipped.py`. Remove the eight moved methods and the
`clipped_reads_dict`/`clipped_reads_bwrd_dict`/`_index` attributes from `ISClipped.__init__`
once moved (they no longer outlive a single `search()` call).

### 3. Out of scope

- Splitting `direction` into separate forward/backward functions or adapters. Deferred (noted
  as architecture review candidate 03 — "Worth exploring," sequenced *after* this ticket
  because splitting is easier once the search logic already lives in its own module).
- Removing the BLAST+ dependency in favor of parsing SAM tags directly. Noted during grilling
  as a future idea, not part of this extraction — the deletion test for it is entirely
  different (BLAST would stop being a dependency at all, not just an injectable one).
- Any change to `search`'s numeric/output behavior relative to today's combination of
  `collect_clipped_reads`/`runblast`/`parseblast`. Byte-for-byte identical, verified below.

## Verification

**Characterization test first** — write `tests/test_clipped_read_search.py` against *today's*
`ISClipped._crtable_ungapped`/`_cl_read_cov_overlap`/`_clboundaries` using
`tests/fake_alignment.py`, pinning current output as a golden value, before doing the move.
Stub-based, consistent with `tests/test_check_blast_ref.py`'s existing convention of not
requiring BLAST+ installed:

- Pure pysam parsing (`_clboundaries`, `_crtable_ungapped`, `_write_cl_fasta`,
  `_cl_read_cov_overlap`) characterized directly against `fake_alignment.py` fixtures — no
  BLAST involved.
- The BLAST round-trip (`parseblast`'s filtering/coordinate-assignment logic) characterized by
  feeding a canned BLAST output table as a fixture, same pattern as
  `tests/test_find_pair.py`/`tests/test_estimate_frequencies.py` (pin against literal
  DataFrames). The new `run_blast` seam is exercised by a fake in the test, never the real
  `blastn` subprocess.

After the move, confirm the same test passes unchanged against `clipped_read_search.search`.

**Manual real-sample diff** — same as tickets 07/09's validation: run the full pipeline on a
real sample in both average and precise mode before and after this change, diff
`ijump_junctions.txt` and `ijump_junction_pairs.txt` byte-for-byte. Must match exactly.

## Done when

- `clipped_read_search.py` exists with a public `search(...)` function and `SearchResult`
  dataclass; `ISClipped` no longer defines the eight moved methods or the
  `clipped_reads_dict`/`clipped_reads_bwrd_dict`/`_index` attributes.
- `ISClipped.run()` calls `clipped_read_search.search` once per direction and assigns results
  onto `self.*`.
- `tests/test_clipped_read_search.py` exists, passes, requires no BLAST+ installation to run.
- Manual real-sample before/after diff of `ijump_junctions.txt` and `ijump_junction_pairs.txt`
  is clean in both modes.
- `pytest` passes from a clean clone.

## Comments
