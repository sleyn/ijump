# 09 — Extract frequency estimation into a standalone `frequency_estimation` module

Status: ready-for-agent
Blocked by: —

## Why

`ISClipped.assess_isel_freq` (`isclipped.py:923-1069`) is the precise pipeline's second
algorithmic core — the frequency-estimation counterpart to the junction-pairing core ticket
08 already moved out. It is ~146 lines behind the trivial interface `assess_isel_freq(self)`,
but the implementation reaches into five separate pieces of `ISClipped` instance state
(`self.pairs_df`, `self.clipped_reads_bwrd`, `self.unclipped_depth`,
`self.cl_read_cov_overlap`, `self.match_lengths`/`read_lengths`/`n_reads_analyzed`). Nothing
can exercise it without first constructing and running the whole pipeline up to that point —
the same gap `.scratch/isclipped-refactor/notes/no-test-safety-net.md` flagged for
`_find_pair` before ticket 06/08 closed it.

Architecture review candidate 01 (`/improve-codebase-architecture`, 2026-08-10). Vocabulary:
`/codebase-design`.

Two supporting facts found during the review, both load-bearing for the scope below:

1. `self.min_match` and `self.av_read_len`, as set inside `assess_isel_freq`
   (`isclipped.py:1000-1001`), are read only later *within that same method*
   (`isclipped.py:1004-1018`). `report_average` (`isclipped.py:1238-1239`) computes the
   identical formulas independently rather than reading what `assess_isel_freq` set. So
   despite living on `self`, these two attributes are effectively local-only inside
   `assess_isel_freq` — nothing outside it depends on the assignment.
2. `ISClipped.fisher_test_clr_number` (`isclipped.py:895-921`) reads `self.min_match`/
   `self.av_read_len` but has zero callers anywhere in the codebase (confirmed by grep,
   including `Test/` and `tests/`) — dead code.

## Scope

### 1. Create the module

Create `frequency_estimation.py` at the repo root (alongside `junction_pairing.py`,
`gff.py`) with a public function:

```python
def estimate_frequencies(pairs_df, clipped_reads_bwrd, unclipped_depth, cl_read_cov_overlap,
                          match_lengths, read_lengths, n_reads_analyzed) -> pd.DataFrame:
    ...
```

Move the body of `assess_isel_freq` (`isclipped.py:923-1069`) into it, adapted from
mutating `self.pairs_df` in place to building and returning a **new** enriched DataFrame
(matching the pattern `junction_pairing.find_pairs` already established — a fresh
DataFrame out, no hidden mutation). `min_match` and `av_read_len` become local variables
inside `estimate_frequencies`, computed exactly as today
(`min(match_lengths)`, `read_lengths / n_reads_analyzed`) — do **not** thread them back out
as return values or as attributes assigned on any object; nothing outside the function
needs them (see "Why," fact 1).

Move the four helpers as module-private functions, no longer methods on `ISClipped`:

| Today | After |
|---|---|
| `ISClipped._calc_freq_precise` (`isclipped.py:794-801`, `@staticmethod`) | `_calc_freq_precise` in `frequency_estimation.py` |
| `ISClipped._read_count_mtx` (`isclipped.py:804-842`, `@staticmethod`) | `_read_count_mtx` in `frequency_estimation.py` |
| `ISClipped._resore_orig_counts` (`isclipped.py:846-854`, `@staticmethod`) | `_resore_orig_counts` in `frequency_estimation.py` |
| `ISClipped._add_total_depth` (`isclipped.py:881-887`, takes unused `self`) | `_add_total_depth` in `frequency_estimation.py`, drop the unused `self` param |

Preserve each helper's body verbatim — relocation, not a rewrite. Keep existing
docstrings/comments.

### 2. Update the call site

`ISClipped.assess_isel_freq` becomes a thin wrapper at its one call site
(`isclipped.py:1176`, inside `run()`):

```python
self.pairs_df = frequency_estimation.estimate_frequencies(
    self.pairs_df, self.clipped_reads_bwrd, self.unclipped_depth,
    self.cl_read_cov_overlap, self.match_lengths, self.read_lengths, self.n_reads_analyzed,
)
```

Add `import frequency_estimation` to `isclipped.py`. Remove `assess_isel_freq` and the four
helper methods from `ISClipped` once moved.

### 3. Delete dead code

Delete `ISClipped.fisher_test_clr_number` (`isclipped.py:895-921`). It has no callers. Once
`min_match`/`av_read_len` stop being set on `self` (per this ticket), it would silently
degrade from "reads a value nothing sets" to "reads a value nothing sets, but now always the
constructor default (150/150)" — a latent trap for whoever revives it later. Deleting it now
avoids stranding it in a more misleading state than it's already in.

## Out of scope

- `count_depth_unclipped` (`isclipped.py:856-880`) stays in `ISClipped`. It does real I/O
  (`self.aln.fetch`) to populate `unclipped_depth` — a data-collection step, not a pure
  computation, and pulling it in would drag `pysam` into what should stay a plain-data
  module.
- Deduplicating the `min_match`/`av_read_len` computation against `report_average`
  (`isclipped.py:1238-1239`). Left exactly as it is today — this ticket only removes the
  redundant assignment onto `self` inside `assess_isel_freq`, it does not touch
  `report_average`'s independent computation. A follow-up ticket if anyone wants to unify
  them.
- Any change to `assess_isel_freq`'s numeric behavior. Byte-for-byte identical output,
  verified per below.

## Verification

**New unit test** — `tests/test_estimate_frequencies.py`. No pre-existing captured-real-data
fixture exists for this function (unlike ticket 06's `_find_pair`, which found one already
pasted into the stale `Test/test_find_pair.py`), and the tiny end-to-end fixture in
`tests/fixtures/` never reaches this code path (it exits at the no-clipped-reads path, per
ticket 03's note). Hand-construct a small `pairs_df` (2-3 rows spanning at least one
orphan-left and one orphan-right row, matching `Position_l`/`Position_r`/`IS_name`/`Chrom`
shapes from `junction_pairing.find_pairs`'s output) plus matching `clipped_reads_bwrd`,
`unclipped_depth`, `cl_read_cov_overlap`, `match_lengths`, `read_lengths`,
`n_reads_analyzed`. Run today's `assess_isel_freq` against it before the move, pin the
output as a golden value. State explicitly in the test docstring that the fixture is
synthetic (hand-built), not observed real data — characterization, not specification.

**Manual real-sample diff** — same as ticket 07's validation: run the full pipeline on a
real sample in precise mode before and after this change, diff `ijump_junction_pairs.txt`
byte-for-byte. Must match exactly.

## Done when

- `frequency_estimation.py` exists with a public `estimate_frequencies(...)` function;
  `ISClipped` no longer defines `assess_isel_freq`'s body, the four helper statics, or
  `fisher_test_clr_number`.
- `ISClipped.run()` calls `frequency_estimation.estimate_frequencies` and assigns the result
  to `self.pairs_df`.
- `tests/test_estimate_frequencies.py` exists, passes, and its docstring states the fixture
  is synthetic.
- Manual real-sample before/after diff of `ijump_junction_pairs.txt` is clean.
- `pytest` passes from a clean clone.

## Comments

**Correction (2026-08-17, review-followups ticket 02):** This ticket's "Done when" required `pytest` to pass from a clean clone. That criterion was not actually met when the ticket was closed: this extraction moved `_read_count_mtx` (and the other helper statics) off `ISClipped` to module level in `frequency_estimation.py` without leaving a delegating alias, which broke `tests/test_no_results_paths.py::test_read_count_mtx_rejects_invalid_orientation` (it called `ISClipped._read_count_mtx(...)`, which no longer existed, so the test raised `AttributeError` instead of the `ValueError` it asserted). Several later tickets recorded that failure as "pre-existing" — it was not; it originated here. Fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.
