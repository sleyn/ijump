# 07 — Collapse ijump.py's orchestration into `ISClipped.run(mode)`

Status: ready-for-agent
Blocked by: —

## Why

`ijump.py` drives `ISClipped` through 16 sequential method calls in a fixed, undocumented
order (`ijump.py:210-331`), each mutating shared instance state that later calls implicitly
depend on (`self.blastout_filtered`, `self.pairs_df`, `self.min_match`, `self.av_read_len`,
`self.boundaries`, `self.clipped_reads`/`clipped_reads_dict`,
`self.cl_read_cov_overlap`/`unclipped_depth`). The interface `ISClipped` exposes today is as
wide as its implementation: the caller has to know the right order, and reaches past the seam
to do it:

- `ijump.py:221,263` builds `self.clipped_reads`/`self.clipped_reads_bwrd` from
  `..._dict` itself — an internal accumulation detail the caller shouldn't need to know.
- `ijump.py:283` overwrites `self.pairs_df` with the output of `filter_pairs`, a free
  function that lives in `ijump.py`, not `isclipped.py`.
- `ijump.py:258` recomputes `av_read_len` inline
  (`is_processing.read_lengths / is_processing.n_reads_analyzed`), duplicating the same
  calculation `ISClipped` does internally later (`isclipped.py:1079`, `1198`).
- `write_empty_outputs` (`ijump.py:34-41`), called from the `NoInsertionsFound` handler,
  already reaches into `ISClipped`'s internal table shapes
  (`ISClipped._cltbl_init()`, `is_processing.sum_by_reg_tbl_init()`, etc.) from outside the
  class — the seam is half-crossed already.

See `.scratch/isclipped-refactor/notes/no-test-safety-net.md` for the broader finding this
was pulled from, and the architecture review that named this candidate 1 (deepen `ISClipped`'s
orchestration interface).

## Scope

Add `ISClipped.run(mode) -> RunResult` (or equivalent status object with at least an
`insertions_found: bool`) as the one entry point `ijump.py` calls after construction:

- Internally sequences the *same* existing methods, in the *same* order, for both
  `average` and `precise` mode — this is a mechanical wrap, not a rewrite. Instance
  attributes (`self.pairs_df`, `self.blastout_filtered`, etc.) stay exactly as they are
  today; do not convert them to local variables in this ticket.
- Absorbs every CSV write currently interleaved in `ijump.py`'s try block: `reads.txt`,
  `ijump_sum_by_reg.txt`, `ijump_report_by_is_reg.txt`, `ijump_junction_pairs.txt`, and
  `reference_regions.tsv` (`ijump.py:224,246,250-251,255,323`).
- Catches `NoInsertionsFound` internally wherever it's currently raised, writes the empty
  outputs itself (folding in today's `write_empty_outputs`), and returns a status instead of
  letting the exception cross the seam. `ijump.py` should end up with no `try/except
  NoInsertionsFound` at all — just logging and `sys.exit(0)` based on the returned status.
- Moves `filter_pairs` (currently a free function in `ijump.py:111-120`) into `isclipped.py`
  — either as a method or a module-level helper `isclipped.py` calls internally — so
  `self.pairs_df` is never assigned from outside the class.
- Moves the `clipped_reads_dict`/`clipped_reads_bwrd_dict` → DataFrame conversions
  (`ijump.py:221,263`) inside `collect_clipped_reads`/`crtable_bwds` respectively.
- Deletes the duplicate `av_read_len` calculation at `ijump.py:258`; `ISClipped` computes it
  once, internally, at the point it's already computed today.
- `ijump.py:main()` shrinks to: build the alignment/args, construct `ISClipped`, call
  `is_processing.run(args.estimation_mode)`, log/exit based on the result, and (unchanged)
  call `create_circos_files()` afterward if `--circos` and mode is average.

## Out of scope

- Converting `self.*` instance attributes to local variables/explicit parameters between the
  now-internal steps — real coupling reduction, but higher risk with no test net beyond
  ticket 06's `_find_pair` characterization test. Leave for a later ticket.
- Any change to the average/precise pipelines' actual behavior, output format, or file
  contents.
- Extracting `_find_pair`/pairing logic — that's ticket 08.

## Validation

No new automated test is required for this ticket. Instead:

- Run `ijump.py` against a real sample in both `average` and `precise` mode before the
  change; save all output files (`reads.txt`, `ijump_sum_by_reg.txt`,
  `ijump_report_by_is_reg.txt`, `ijump_junction_pairs.txt`, `ijump.log` structure aside) and
  `--circos` output.
- Run again after the change with the same inputs; diff every output file byte-for-byte
  against the pre-change run. They must match exactly.
- Also run once in a scenario that hits `NoInsertionsFound` (or verify the code path by
  inspection if no such sample is on hand) to confirm the empty-output files still match.

## Done when

- `ijump.py:main()` calls `is_processing.run(args.estimation_mode)` and contains no direct
  reference to `ISClipped` internal attributes or `filter_pairs`.
- `filter_pairs` lives in `isclipped.py`, not `ijump.py`.
- `write_empty_outputs` and the `NoInsertionsFound` try/except in `ijump.py` are gone;
  `run()` handles both internally.
- Before/after output diff (manual run, both modes) is clean.
- `pytest` still passes (existing suite, including ticket 06's `_find_pair` test).

## Comments
