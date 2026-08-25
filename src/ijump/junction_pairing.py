# Junction pairing: pair left and right insertion junctions into IS element positions.

import numpy as np
import pandas as pd


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

    # Check if both left and right junctions present. If not - process just present
    # part of junctions.
    if pos_l.size == 0 or pos_r.size == 0:
        n_pairs = pos_l.size + pos_r.size
        pairs_df = pd.DataFrame(
            {
                "Position_l": [0] * n_pairs,
                "Position_r": [0] * n_pairs,
                "Count_mapped_to_IS_l": [0] * n_pairs,
                "Count_mapped_to_IS_r": [0] * n_pairs,
                "Chrom": chrom,
            }
        )

        if pos_r.size == 0:
            for pos_l_index, pos in enumerate(pos_l):
                pairs_df.iloc[pos_l_index, :] = [pos, 0, pos_l_count[pos_l_index], 0, chrom]
        else:
            for pos_r_index, pos in enumerate(pos_r):
                pairs_df.iloc[pos_r_index, :] = [pos, 0, pos_r_count[pos_r_index], 0, chrom]

        return pairs_df

    # Store close positions in the matrix where rows are left positions and columns
    # are right positions.
    # The value is 1 if two positions are closer then max_is_dup_len value.
    closeness_matrix = np.zeros((pos_l.size, pos_r.size))

    # Check if any position close to the contig ends.
    if pos_r[-1] - pos_l[0] > chrom_len / 2:
        if chrom_len - (pos_r[-1] - pos_l[0]) <= max_is_dup_len:
            closeness_matrix[0, -1] = 1

    if pos_l[-1] - pos_r[0] > chrom_len / 2:
        if chrom_len - (pos_l[-1] - pos_r[0]) <= max_is_dup_len:
            closeness_matrix[-1, 0] = 1

    # Populate closeness matrix.
    for pos_index, pos in enumerate(pos_l):
        closeness_matrix[pos_index] = (
            np.ones_like(pos_r) * (np.abs(pos_r - pos) < max_is_dup_len)
        ).astype(np.intp)

    # Assign clusters and sort in each cluster by junction representation in descending order.

    # Build dataframe to populate pairs.
    # We will use maximum number of rows (if all positions do not have pairs).
    n_pairs = np.sum(closeness_matrix.shape)

    pairs_df = pd.DataFrame(
        {
            "Position_l": [0] * n_pairs,
            "Position_r": [0] * n_pairs,
            "Count_mapped_to_IS_l": [0] * n_pairs,
            "Count_mapped_to_IS_r": [0] * n_pairs,
            "Chrom": chrom,
        }
    )

    # Build clusters of close positions.
    # Clusters are attributed to the left joints.
    cluster_ids = np.zeros(len(closeness_matrix))
    cluster_cur_id = 0
    column_index = 0
    # Itrerate through all closeness_matrix columns or until all left joints will be
    # assigned to clusters.
    while not (column_index >= closeness_matrix.shape[1] or np.all(cluster_ids > 0)):
        # Check if the column not zero (orphan right position).
        if closeness_matrix[:, column_index].any():
            # If any left position has several right positions in proximity
            # unite clusters.
            if np.any(closeness_matrix[:, column_index][cluster_ids > 0] == 1):
                cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id
            else:
                # If cluster is first or clusters do not overlap add cluster id.
                cluster_cur_id += 1
                cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id

        column_index += 1

    # Sort each cluster.
    for cluster_id in np.unique(cluster_ids[cluster_ids > 0]):
        # Sort positions sub-list.
        cluster_pos_l = pos_l[np.where(cluster_ids == cluster_id)]
        cluster_pos_l = cluster_pos_l[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        pos_l[np.where(cluster_ids == cluster_id)] = cluster_pos_l

        # Sort closeness sub-matrix.
        cluster_closeness_matrix = closeness_matrix[np.where(cluster_ids == cluster_id)]
        cluster_closeness_matrix = cluster_closeness_matrix[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        closeness_matrix[np.where(cluster_ids == cluster_id), :] = cluster_closeness_matrix

        # Sort counsts sub-list.
        cluster_pos_l_count = pos_l_count[np.where(cluster_ids == cluster_id)]
        cluster_pos_l_count = cluster_pos_l_count[
            np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]
        ]
        pos_l_count[np.where(cluster_ids == cluster_id)] = cluster_pos_l_count

    # Collect right indexes.
    pos_r_orphan = np.arange(pos_r.size)

    # Populate pairs table.
    for pos_l_index, pos_l_cur in enumerate(pos_l):
        if np.sum(closeness_matrix[pos_l_index, :]):
            # Find right index that has minimum difference in counts.
            # Penalize non-cluster items difference by 10000.
            pos_r_index = np.argmin(
                np.abs(pos_r_count - pos_l_count[pos_l_index])
                + ~(closeness_matrix[pos_l_index, :] == 1) * 10000
            )
            pairs_df.iloc[pos_l_index, :] = [
                pos_l_cur,
                pos_r[pos_r_index],
                pos_l_count[pos_l_index],
                pos_r_count[pos_r_index],
                chrom,
            ]

            closeness_matrix[:, pos_r_index] = 0

            pos_r_orphan[pos_r_index] = -1

        # Write orhphan peaks.
        else:
            pairs_df.iloc[pos_l_index, :] = [pos_l_cur, 0, pos_l_count[pos_l_index], 0, chrom]

    df_offset = len(pos_l)

    # Add right orphan peaks.
    for shift, pos_r_index_orphan in enumerate(pos_r_orphan[pos_r_orphan != -1]):
        pairs_df.iloc[df_offset + shift, :] = [
            0,
            pos_r[pos_r_index_orphan],
            0,
            pos_r_count[pos_r_index_orphan],
            chrom,
        ]

    # Remove empty rows.
    pairs_df = pairs_df.query("Position_l > 0 or Position_r > 0")

    return pairs_df
