# 06 — Characterization test for `_find_pair`

Status: ready-for-agent
Blocked by: —

## Why

`ISClipped._find_pair` (`isclipped.py:667-...`, a `@staticmethod`) is the junction-pairing
core — ~155 lines, the part of `ISClipped` most likely to change shape under any future
refactor and the hardest to reason about by eye. It currently has zero executable
verification.

The only existing test touching it, `Test/test_find_pair.py`, is not a real test:

1. **It tests a copy, not the code.** The file pastes its own definition of `_find_pair` at
   line 5 rather than importing `ISClipped._find_pair` from `isclipped.py`. The shipped
   implementation is never loaded. The copy has already drifted from the real one — the
   pasted version builds its output dict with a `Count_l`/`Count_r` column pair, while the
   shipped version at `isclipped.py:672-675` uses `Count_mapped_to_IS_l`/`_r`.
2. **It has no assertions.** It calls the function, binds the result to `x`, and prints
   `'OK'` unconditionally. It exits 0 whatever the function returns.

See `.scratch/isclipped-refactor/notes/no-test-safety-net.md` for the original finding.

## Scope

Add a tracked test (in `tests/`, alongside the suite from ticket 01) that:

- Imports the real `ISClipped._find_pair` from `isclipped.py` — no pasted copy.
- Calls it with the real-data fixture already sitting in `Test/test_find_pair.py:165-171`
  (`pos_l`, `pos_r`, `pos_l_count`, `pos_r_count`, `chrom_len = 3909467`,
  `max_is_dup_len = 20`, and the contig name) — captured from actual data, so it exercises
  real coordinate ranges.
- Pins the returned `pairs_df` as a golden value and asserts the result equals it (e.g.
  `pandas.testing.assert_frame_equal`), including column names — the `Count_mapped_to_IS_l`
  discrepancy above should be caught by column-name assertions, not just values.

This is characterization, not specification: the goal is "did I change behaviour", not "is
this correct". Do not attempt to fix or validate the algorithm's correctness here.

Once the new test is in place and passing, delete `Test/test_find_pair.py` — its pasted copy
is now redundant and misleading (it doesn't exercise the shipped code).

## Constraints

- No BAM, BLAST, or GFF input needed — `_find_pair` takes plain numpy arrays and a static
  method call needs no `ISClipped` instance.
- Keep the fixture inline in the test file (it's a handful of small arrays); no new fixture
  files needed under `tests/fixtures/`.

## Out of scope

- End-to-end snapshot testing over the full `Test/` sample data (`Sample.bam`, the
  A. baumannii references) — the note's broader suggestion, but a separate, larger effort
  given `Test/` is gitignored and holds ~8.7 GB of data that can't be committed as-is.
- Any change to `_find_pair`'s behavior or its callers.

## Done when

- `tests/` has a test that imports and calls the real `ISClipped._find_pair`.
- The test asserts a pinned golden output (values and column names), not just "doesn't
  crash".
- `pytest` passes from a clean clone with no bioinformatics tooling installed.
- `Test/test_find_pair.py` is deleted.

## Comments

Implemented: added `tests/test_find_pair.py`, importing `ISClipped._find_pair` directly
from `isclipped.py` and calling it with the fixture from the old
`Test/test_find_pair.py:165-171`. The golden `pairs_df` was captured by running the real
implementation once and pinned with `pandas.testing.assert_frame_equal` (columns and
values), so it does catch the `Count_mapped_to_IS_l`/`Count_l` naming discrepancy noted
above. `Test/test_find_pair.py` deleted.

Note: this shell's Python environment lacks `pysamstats`, so `isclipped.py` (and therefore
every test importing it, including the ticket 01/02 suite) fails to collect here — a
pre-existing environment gap, not introduced by this change. Verified the new test's logic
by temporarily stubbing `pysamstats`/`Bio.Blast.Applications`/`sklearn.cluster` in
`sys.modules` for a one-off run against the real `_find_pair`; the stub was not committed.
