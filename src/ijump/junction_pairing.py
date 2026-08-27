# Junction pairing: pair left and right insertion junctions into IS element positions.

import numpy as np
import pandas as pd

# What a row carries on the side that has no junction.
#
# Positions here are 0-based, so 0 is the first base of a contig -- a real
# coordinate, and a reachable one: an origin-spanning element leaves a locus at
# the very start of its contig, and a read clipped there reports its junction at
# 0. Absence therefore cannot be spelled 0 without swallowing that junction.
#
# The written file is unaffected. It is 1-based, where 0 is not a position at
# all, so absence is written there as 0 -- unambiguously -- by
# ``isclipped.convert_zero_one_base``.
NO_JUNCTION = -1

# What `pos_r_orphan` below holds for a right junction that found a partner. It
# is an index into pos_r, not a coordinate -- the same literal as NO_JUNCTION and
# a wholly unrelated meaning, so it is named rather than left to read as one.
PAIRED = -1


# A pairs-table row is written as a bare positional list, so which column a value
# lands in is carried by nothing but its place in that list. These three builders
# name the columns once so no write site has to get the order right.
def _pair_row(pos_l, pos_r, count_l, count_r, chrom):
    return [pos_l, pos_r, count_l, count_r, chrom]


def _left_orphan_row(pos, count, chrom):
    """A left junction with no right partner."""
    return _pair_row(pos, NO_JUNCTION, count, 0, chrom)


def _right_orphan_row(pos, count, chrom):
    """A right junction with no left partner."""
    return _pair_row(NO_JUNCTION, pos, 0, count, chrom)


def _empty_pairs_frame(n_rows, chrom):
    """A frame of ``n_rows`` rows for the writers below to fill in.

    Positions start at NO_JUNCTION so that a row nobody writes to is absent on
    both sides, which is what lets the main path tell its unused rows from a
    junction at position 0.
    """
    return pd.DataFrame(
        {
            "Position_l": [NO_JUNCTION] * n_rows,
            "Position_r": [NO_JUNCTION] * n_rows,
            "Count_mapped_to_IS_l": [0] * n_rows,
            "Count_mapped_to_IS_r": [0] * n_rows,
            "Chrom": chrom,
        }
    )


# Make clusters of left and right insertions junctions from positions.
# Outputs table of right and left positions of of junctions pairs with counts of clipped
# reads accounted to each junction.
# Similar results could be achieved by KNN search, but this algorithm shows slightly
# better performance on tests.
def find_pairs(pos_l, pos_r, pos_l_count, pos_r_count, chrom_len, max_is_dup_len, chrom):
    # The cluster sort below writes back into pos_l and pos_l_count, so work on
    # copies: the caller passes Series.to_numpy() results, which pandas 3 hands
    # out non-writeable (pandas 2 returned writeable copies, which is why this
    # only surfaced on the upgrade). Copying also keeps the caller's arrays out
    # of the reordering, which was never intentional. pos_r and pos_r_count are
    # only read today, but are copied too so the whole signature carries one
    # rule -- inputs are never touched -- rather than a per-argument exception.
    pos_l = np.array(pos_l, copy=True)
    pos_r = np.array(pos_r, copy=True)
    pos_l_count = np.array(pos_l_count, copy=True)
    pos_r_count = np.array(pos_r_count, copy=True)

    if pos_l.size == 0 or pos_r.size == 0:
        n_pairs = pos_l.size + pos_r.size
        pairs_df = _empty_pairs_frame(n_pairs, chrom)

        # Nothing to pair, so every junction is an orphan on its own side --
        # written the same way the main path below writes the orphans it finds.
        if pos_r.size == 0:
            for pos_l_index, pos in enumerate(pos_l):
                pairs_df.iloc[pos_l_index, :] = _left_orphan_row(
                    pos, pos_l_count[pos_l_index], chrom
                )
        else:
            for pos_r_index, pos in enumerate(pos_r):
                pairs_df.iloc[pos_r_index, :] = _right_orphan_row(
                    pos, pos_r_count[pos_r_index], chrom
                )

        return pairs_df

    # Rows are left positions, columns are right positions. The value is 1 if two
    # positions are closer than max_is_dup_len, where "closer" accounts for
    # wraparound at the contig's circular boundary: a left junction near position
    # 0 and a right junction near chrom_len can be a genuine close pair (an
    # origin-spanning insertion), so distance is measured both ways round the
    # contig and the shorter one is used.
    closeness_matrix = np.zeros((pos_l.size, pos_r.size))

    for pos_index, pos in enumerate(pos_l):
        linear_dist = np.abs(pos_r - pos)
        circular_dist = np.minimum(linear_dist, chrom_len - linear_dist)
        closeness_matrix[pos_index] = (circular_dist < max_is_dup_len).astype(np.intp)

    # Worst case every junction on both sides is an orphan, so size for that.
    n_pairs = np.sum(closeness_matrix.shape)

    pairs_df = _empty_pairs_frame(n_pairs, chrom)

    # Clusters of close positions, attributed to the left joints.
    cluster_ids = np.zeros(len(closeness_matrix))
    cluster_cur_id = 0
    column_index = 0
    while not (column_index >= closeness_matrix.shape[1] or np.all(cluster_ids > 0)):
        if closeness_matrix[:, column_index].any():
            # A left position already in a cluster that is also close to this
            # right position merges into it; otherwise this right position starts
            # a new cluster.
            if np.any(closeness_matrix[:, column_index][cluster_ids > 0] == 1):
                cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id
            else:
                cluster_cur_id += 1
                cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id

        column_index += 1

    # Sort each cluster's left positions (and the count/closeness rows that go
    # with them) by descending clipped-read count, highest-confidence first.
    for cluster_id in np.unique(cluster_ids[cluster_ids > 0]):
        cluster_pos_l = pos_l[np.where(cluster_ids == cluster_id)]
        cluster_pos_l = cluster_pos_l[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        pos_l[np.where(cluster_ids == cluster_id)] = cluster_pos_l

        cluster_closeness_matrix = closeness_matrix[np.where(cluster_ids == cluster_id)]
        cluster_closeness_matrix = cluster_closeness_matrix[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        closeness_matrix[np.where(cluster_ids == cluster_id), :] = cluster_closeness_matrix

        cluster_pos_l_count = pos_l_count[np.where(cluster_ids == cluster_id)]
        cluster_pos_l_count = cluster_pos_l_count[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        pos_l_count[np.where(cluster_ids == cluster_id)] = cluster_pos_l_count

    pos_r_orphan = np.arange(pos_r.size)

    for pos_l_index, pos_l_cur in enumerate(pos_l):
        if np.sum(closeness_matrix[pos_l_index, :]):
            # Find right index that has minimum difference in counts.
            # Penalize non-cluster items difference by 10000.
            pos_r_index = np.argmin(
                np.abs(pos_r_count - pos_l_count[pos_l_index])
                + ~(closeness_matrix[pos_l_index, :] == 1) * 10000
            )
            pairs_df.iloc[pos_l_index, :] = _pair_row(
                pos_l_cur,
                pos_r[pos_r_index],
                pos_l_count[pos_l_index],
                pos_r_count[pos_r_index],
                chrom,
            )

            closeness_matrix[:, pos_r_index] = 0

            pos_r_orphan[pos_r_index] = PAIRED

        else:
            pairs_df.iloc[pos_l_index, :] = _left_orphan_row(
                pos_l_cur, pos_l_count[pos_l_index], chrom
            )

    df_offset = len(pos_l)

    for shift, pos_r_index_orphan in enumerate(pos_r_orphan[pos_r_orphan != PAIRED]):
        pairs_df.iloc[df_offset + shift, :] = _right_orphan_row(
            pos_r[pos_r_index_orphan], pos_r_count[pos_r_index_orphan], chrom
        )

    # Drop the rows nobody wrote to. The frame is sized for the worst case --
    # every junction on both sides an orphan -- and the writers above fill it
    # from the front: one row per left junction, then one per right orphan. So
    # the used rows are the leading `df_offset + n_right_orphans` of them, and
    # slicing takes exactly those. Selecting them by value instead ("keep rows
    # with a position above zero") is what dropped a junction at position 0,
    # which is a coordinate and not an absence.
    n_written = df_offset + int(np.sum(pos_r_orphan != PAIRED))

    return pairs_df.iloc[:n_written].copy()
