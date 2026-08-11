# 08 — Move `_find_pair` into a standalone `junction_pairing` module

Status: ready-for-agent
Blocked by: —

## Why

`ISClipped._find_pair` (`isclipped.py:667-819`, a `@staticmethod`) is the junction-pairing
core — ~155 lines, the hardest algorithm in the file — but it is already a fully
parameter-only function: `pos_l, pos_r, pos_l_count, pos_r_count, chrom_len, max_is_dup_len,
chrom`. It reads and writes no `self`. Ticket 06 already gave it a passing characterization
test (`tests/test_find_pair.py`) that imports and calls it directly.

The deletion test and interface-shape work for this candidate are effectively done — what's
left is placement. Reaching this algorithm today still requires importing all 1332 lines of
`isclipped.py`, and its `_find_pair` name and location signal "private implementation detail
of `ISClipped`" when it's actually an independent, reusable piece of the pairing logic.
`CONTEXT.md`'s Precise mode entry already names this domain concept: "attempts to resolve
each insertion event separately by pairing left and right junctions."

## Scope

- Create `junction_pairing.py` at the repo root (alongside `isclipped.py`, `gff.py`, etc.)
  with a public function:

  ```python
  def find_pairs(pos_l, pos_r, pos_l_count, pos_r_count, chrom_len, max_is_dup_len, chrom):
      ...
  ```

  Move the body of `ISClipped._find_pair` (`isclipped.py:667-819`) into it verbatim —
  relocation, not a rewrite. Keep the existing docstring/comments.
- Remove `_find_pair` from `ISClipped`; update its one call site,
  `search_insert_pos` (`isclipped.py:854`), to call `find_pairs(...)` from the new module
  instead of `self._find_pair(...)`. Add `from junction_pairing import find_pairs` (or
  `import junction_pairing`) to `isclipped.py`.
- Update `tests/test_find_pair.py` to import `find_pairs` from `junction_pairing` instead of
  `ISClipped` from `isclipped`. The golden fixture and assertion stay unchanged — same
  inputs, same expected `pairs_df`.

## Out of scope

- Any change to the algorithm's behavior, including the `pos_l`/`pos_l_count` in-place
  mutation it currently does internally (`isclipped.py:757,763,769` equivalents) — preserve
  exactly.
- Renaming or restructuring `search_insert_pos`'s surrounding loop (`isclipped.py:822-869`)
  beyond the one call-site change.
- Ticket 07's orchestration work.

## Done when

- `junction_pairing.py` exists with a public `find_pairs(...)` function; `isclipped.py` no
  longer defines `_find_pair`.
- `search_insert_pos` calls `junction_pairing.find_pairs`.
- `tests/test_find_pair.py` imports from `junction_pairing`, not `isclipped`, and still
  passes with the same golden output.
- `pytest` passes from a clean clone (module import structure aside — see ticket 06's
  comment on the `pysamstats` environment gap, not something this ticket needs to fix).

## Comments
