# 03 — Fix the 'presice' typo: IS pos is never converted to 1-base

Status: ready-for-agent
Blocked by: 01

## The bug

`ijump.py:39`:

```python
if est_mode == 'presice':
    junc_tbl_copy['IS pos'] += 1
```

The CLI value is `precise` (ijump.py:113, checked at ijump.py:220). `'presice'` is a typo,
so this branch has **never executed**. The `IS pos` column in `ijump_junctions.txt` has
always been left 0-based while `Position` on the next line is converted to 1-based
(ijump.py:41).

Every `ijump_junctions.txt` ever produced in precise mode has a `IS pos` column off by one
relative to `Position`.

## Scope

Fix the comparison to `'precise'`.

Then decide — and record the decision in the ticket comments — whether `IS pos` should be
1-based at all. The typo means the intended behaviour was never observed, so "what it was
supposed to do" is an assumption, not an established fact. Check against `docs/Precise.md`
and the column's use in `combine_results.py` before committing to the conversion.

Prefer comparing against a shared constant or enum over a bare string literal, so the next
typo is a `NameError` rather than a silently dead branch. Note ijump.py:201 and :220 compare
the same values.

## Why it is a separate ticket

Deliberately kept out of ticket 02. Fixing it shifts data values inside
`ijump_junctions.txt`, so a before/after diff of ticket 02's refactor would no longer be
empty — making an intended fix indistinguishable from a refactor regression.

## Verification

A regression test asserting `IS pos` and `Position` use the same coordinate base in
precise-mode output. Uses the fixture from ticket 01 — but note that fixture produces no
clipped reads, so it never reaches junction output. Either extend it with a variant that
does produce junctions, or test `check_junctions_presence` directly with a constructed
junctions table (the cheaper option — it is a plain function over a DataFrame,
ijump.py:35).

## User impact

This changes output values. Worth a note in the README or a version bump, since users
comparing new runs against archived ones will see `IS pos` shift by 1.

## Comments

Decision: yes, convert `IS pos` to 1-based in precise mode, matching `Position`.

Both columns are populated from the same coordinate system: `IS pos` comes from
`clipped_reads['clip_position']` or `blastout_filtered['pos_in_ref']` (isclipped.py:632,
:640), and `Position` comes from the other of that same pair (isclipped.py:628, :637) — both
ultimately 0-based positions derived from pysam/BLAST hit coordinates. `docs/Precise.md` and
`docs/Average.md` describe `IS pos` and `Position` as a matched pair of coordinates in the
same junction row, with no indication one should stay 0-based while the other is 1-based.
`combine_results.py` does not read `ijump_junctions.txt` or the `IS pos` column at all (it
consumes `ijump_junction_pairs.txt`'s `Position_l`/`Position_r` instead), so there is no
downstream consumer whose expectations constrain the choice either way.

Implemented: fixed `'presice'` -> `'precise'`, replaced the bare string literals for
`--estimation_mode` values with an `EstimationMode` str enum (`ijump.py`) used at every
comparison site, added a regression test (`tests/test_junctions_coordinate_base.py`) driving
`check_junctions_presence` directly against a constructed junctions table, bumped the
version to 1.0.4, and noted the output-value shift in the README.