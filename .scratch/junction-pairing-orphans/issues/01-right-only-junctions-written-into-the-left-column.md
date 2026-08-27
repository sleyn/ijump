# 01 — Right-only junctions are written into the left column

**What to build:** A cluster whose junctions on a contig are all right-handed is reported
as right orphans carrying their real counts and frequencies, instead of as left
orphans with zero counts and frequency `0.0`.

`find_pairs` returns early when one side has no junctions at all, and that early return has
two arms. The left-only arm writes the position and count into the left pair of columns,
which is right. The right-only arm writes them into the *same* left pair, leaving both
right columns zero — a copy-paste slip, not a convention: the main path's own right-orphan
writer puts a right junction in `Position_r`/`Count_mapped_to_IS_r`, and that is the shape
this arm should produce too.

Downstream, frequency estimation reads a right junction's evidence out of the right
columns, finds zeros there, and reports `0.0`. So the defect does not merely mislabel a
row — it silently zeroes a real detection.

**This is live, not latent.** On the reference sample five (cluster, contig) pairings are
one-sided and two of them are right-only: `ISAba12` and `ISAba53` on `NODE_2`. Those are
the seven rows in the pinned precise-mode golden that report `Frequency 0.0` at positions
451, 145216, 145218, 145219, 145220, 145221 and 145222.

Found while implementing isfinder-annotation 06, which regrouped `IS17`'s two `NODE_1`
junctions under `ISAba12`. That gave the group left junctions as well, so it stopped taking
the early return, took the ordinary path, and the two rows went from `0.0` to real
frequencies. The regrouping routed around the bug for one group; it did not fix it.

Note that `filter_pairs` matches on both position columns, so moving a coordinate from one
to the other can change which rows survive filtering. The golden diff is the place to check
that.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] Right-only input is written as right orphans: position in `Position_r`, count in `Count_mapped_to_IS_r`, both left columns zero
- [x] The left-only arm is unchanged
- [x] Unit coverage pins both arms of the early return — today's `find_pairs` test exercises only the two-sided path
- [x] The precise-mode end-to-end golden is re-pinned, and the diff is confined to the seven `ISAba12`/`ISAba53` rows on `NODE_2`
- [x] Those seven rows carry non-zero frequencies where the underlying evidence is non-zero

## Comments

- The golden diff is the seven predicted rows and nothing else; every other row is
  byte-identical and the row count is unchanged. `filter_pairs` kept all seven, so moving
  the coordinate between columns did not change what survives filtering.
- The zeroed frequencies were not near-zero. `N_clipped_l` was 0.0 for these rows because a
  right junction's clipped reads live in the *right* count matrix, and the left lookup
  found nothing there. Read from the right matrix, position 145221 has 284 clipped reads
  against 2 unclipped, and reports frequency **0.983** where it used to report 0.0 — a
  near-fixed ISAba12 insertion the pipeline was reporting as absent. The other six land
  between 0.0034 and 0.207.
- The precise-mode `ijump_junctions.txt` golden is unchanged; only the pairs table moved.
- Downstream, `combine_results` renames `Position_l`/`Position_r` to `Start`/`Stop` and has
  its own handling for rows where one of the two is 0. The seven rows now reach it in the
  shape it expects for a right orphan instead of masquerading as left orphans.
- The five write sites in `find_pairs` built their rows as bare positional lists, which is
  how the swapped columns went unnoticed. They now go through `_pair_row` /
  `_left_orphan_row` / `_right_orphan_row`, so the column order is written down once.
