# 03 — Wraparound junction pairing only checked the single most-extreme pair

**What to build:** Every left/right junction pair close to each other through
a contig's circular boundary is paired, not just the single smallest-vs-largest
pair.

## Why

`junction_pairing.find_pairs` special-cased circular-boundary proximity (an
origin-spanning insertion's left junction near position 0 and right junction
near `chrom_len`, or vice versa — see `CONTEXT.md`'s **Origin-spanning
element** entry, a first-class supported case) by checking only `pos_l[0]`
vs `pos_r[-1]` and `pos_l[-1]` vs `pos_r[0]`: the smallest-vs-largest position
on each side, since `isclipped.search_insert_pos` always passes these arrays
sorted ascending. Any other pair of junctions close to each other through the
origin was invisible — the ordinary closeness loop measures plain linear
distance (`np.abs(pos_r - pos)`), which is huge for two positions on opposite
sides of the origin that are actually close going the other way round.

Concrete counterexample that exposed it: `chrom_len=1,000,000`,
`max_is_dup_len=20`, left junctions `[1, 5, 999995]`, right junctions
`[3, 999990, 999998]`. Only `(1, 999998)` — the extreme pair — was flagged as
close. `(5, 999990)`, genuinely 15bp apart through the origin, was invisible
to the old code and would have been reported as two orphans (degraded/zeroed
evidence) instead of one paired insertion.

The reference test genome already carries one real origin-spanning locus
(`IS17_2` starting at base 2 of `NODE_2`, established in
`.scratch/junction-pairing-orphans/issues/02-position-zero-means-both-a-coordinate-and-absent.md`'s
Comments), which is why this scenario is a live concern for this repo's own
sample data, not a purely theoretical one.

Found during the same `/grilling` session as tickets 01 and 02, auditing
`find_pairs` for correctness at the user's request.

## Scope

- Replaced both special-cased extreme-pair checks with a single circular
  distance computation — `min(linear_dist, chrom_len - linear_dist) <
  max_is_dup_len` — inside the same vectorized loop that builds the rest of
  the closeness matrix, so every pair is covered rather than just the two
  extremes.
- Normalized the comparison to strict `<` throughout. The two special cases
  previously used `<=`, inconsistent with the main loop's `<`, with nothing
  in the docs or comments explaining an intentional difference.

## Out of scope

- The frequency-formula fixes — see tickets 01 and 02, different functions.
- The clustering/merge logic later in `find_pairs` (the `cluster_ids`
  transitive-merge loop) — noted as a candidate for further audit during the
  session but not investigated to the same depth; not touched here.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] Circular distance, not just the extreme pair, drives wraparound
      closeness.
- [x] `<` used consistently for both linear and circular closeness.
- [x] `tests/test_find_pair.py` gained
      `test_wraparound_pairs_every_close_junction_not_just_the_extreme_pair`,
      pinning the corrected behaviour on the counterexample above (all three
      pairs now correctly matched).
- [x] `tests/goldens/e2e/precise/report/ijump_junction_pairs.txt` re-pinned;
      row count unchanged (44 before and after) — the real reference sample
      doesn't happen to carry multiple near-boundary junctions on one contig,
      so the fix is confirmed behaviour-preserving there (reachability
      checked via the synthetic counterexample, not assumed on faith) while
      being provably correct on the counterexample itself, the same pattern
      `junction-pairing-orphans/02` used.
- [x] `ruff`/`ruff format`/`mypy` clean on `src/ijump/`.

## Comments

Landed alongside tickets 01 and 02 in this same directory, from one audit
session. Full non-e2e suite (190 passed) and full e2e tier (13 passed) green.
The clustering/merge logic noted as out of scope above is a candidate for a
future ticket if someone wants to dig further into `find_pairs`.
