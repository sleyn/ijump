# 01 — Right-only junction groups are written into the left column

**What to build:** A cluster whose junctions on a contig are all right-handed is reported
as right-junction orphans carrying their real counts and frequencies, instead of as left
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

**This is live, not latent.** On the reference sample five (cluster, contig) groups are
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

**Status:** ready-for-agent

- [ ] A right-only group is written as right-junction orphans: position in `Position_r`, count in `Count_mapped_to_IS_r`, both left columns zero
- [ ] The left-only arm is unchanged
- [ ] Unit coverage pins both arms of the early return — today's `find_pairs` test exercises only the two-sided path
- [ ] The precise-mode end-to-end golden is re-pinned, and the diff is confined to the seven `ISAba12`/`ISAba53` rows on `NODE_2`
- [ ] Those seven rows carry non-zero frequencies where the underlying evidence is non-zero
