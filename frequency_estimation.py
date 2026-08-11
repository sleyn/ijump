# Frequency estimation: assess frequency of IS insertions in the population from
# paired junction positions, clipped-read counts, and unclipped/overlap depth.

import numpy as np
import pandas as pd


# Calculate frequency of IS insertion based on frequencies of boundaries junctions.
def _calc_freq_precise(freq_l, freq_r):
    if freq_l == 0:
        return freq_r
    elif freq_r == 0:
        return freq_l
    else:
        return (freq_l + freq_r) / 2


# Make read count matrices.
def _read_count_mtx(pairs_df, orientation):
    if orientation == 'left':
        pos_df = pairs_df.query('Position_l > 0').copy()
        pos_df = pos_df.rename(columns={'Position_l': 'Position', 'Count_mapped_to_IS_l': 'Count'})
    elif orientation == 'right':
        pos_df = pairs_df.query('Position_r > 0').copy()
        pos_df = pos_df.rename(columns={'Position_r': 'Position', 'Count_mapped_to_IS_r': 'Count'})
    else:
        raise ValueError('the parameter should be "left" or "right"')

    # Dictionary to translate positions (Contig name/Coordinate) to matrix row indeces.
    pos = {}
    for chrom in pos_df.Chrom.unique():
        pos[chrom] = {}

    i = 0

    for pos_row in pos_df[['Position', 'Chrom']].drop_duplicates().itertuples():
        pos[pos_row.Chrom][pos_row.Position] = i
        i += 1

    # Dictionary to translate IS names to matrix row indeces.
    is_names = dict(
        zip(
            pos_df['IS_name'].unique().tolist(),
            [i for i in range(pos_df['IS_name'].unique().size)]
        )
    )

    counts = np.zeros((
        len(pos_df[['Position', 'Chrom']].drop_duplicates()),
        len(is_names)
    ))

    for row in pos_df.itertuples():
        counts[pos[row.Chrom][row.Position], is_names[row.IS_name]] = row.Count

    return pos, is_names, counts


# Translate matrix of read count proportions to the original read counts.
# (calculated from the second pass of clipped reads collection)
def _resore_orig_counts(counts_mtx, original_rc_counts, pos_dict):
    counts_mtx /= counts_mtx.sum(1).reshape(-1, 1)

    for chr in pos_dict.keys():
        for junct_pos in pos_dict[chr].keys():
            counts_mtx[pos_dict[chr][junct_pos]] *= original_rc_counts[chr].get(junct_pos, 0)

    return counts_mtx


def _add_total_depth(depth_l, depth_r):
    if depth_l == 0:
        return depth_r
    elif depth_r == 0:
        return depth_l
    else:
        return (depth_l + depth_r) / 2


# Assess frequency of the insertion in population.
def estimate_frequencies(pairs_df, clipped_reads_bwrd, unclipped_depth, cl_read_cov_overlap,
                          match_lengths, read_lengths, n_reads_analyzed) -> pd.DataFrame:
    # Setup calculation of number of reads supporting each position count
    # for each IS element.
    # 1: Count reads for right and left positions that came directly from positions.
    # Caveat - they do not have information about corresponding IS elements.
    pairs_df = pairs_df.copy()

    original_rc_l_df = clipped_reads_bwrd. \
        query('clip_position == "left"'). \
        groupby(['junction_in_read', 'IS_chrom'], as_index=False)['Read name']. \
        count(). \
        rename(columns={'Read name': 'Count', 'IS_chrom': 'Chrom'})

    # Pre-populate one dict entry per chrom seen in pairs_df so that every chrom
    # _resore_orig_counts will look up is guaranteed present, even with zero
    # clipped reads recorded here. Originally pre-populated from ISClipped.ref_len
    # (the full reference chrom set); pairs_df's own chroms are a subset of that
    # and are the only ones ever looked up, so this is behaviorally equivalent
    # without needing ref_len as a parameter.
    original_rc_l = {chrom: {} for chrom in pairs_df['Chrom'].unique()}
    for rc_row_l in original_rc_l_df.itertuples():
        original_rc_l.setdefault(rc_row_l.Chrom, {})[rc_row_l.junction_in_read] = rc_row_l.Count

    original_rc_r_df = clipped_reads_bwrd. \
        query('clip_position == "right"'). \
        groupby(['junction_in_read', 'IS_chrom'], as_index=False)['Read name']. \
        count(). \
        rename(columns={'Read name': 'Count', 'IS_chrom': 'Chrom'})

    original_rc_r = {chrom: {} for chrom in pairs_df['Chrom'].unique()}
    for rc_row_r in original_rc_r_df.itertuples():
        original_rc_r.setdefault(rc_row_r.Chrom, {})[rc_row_r.junction_in_read] = rc_row_r.Count

    # 2: Make matrix for left positions.
    pos_l, is_names_l, counts_l = _read_count_mtx(pairs_df, 'left')

    # 3: Make matrix for right positions.
    pos_r, is_names_r, counts_r = _read_count_mtx(pairs_df, 'right')

    # Calculate proportions of reads for each IS for each conflicting position
    # and split reads supporting position from the clipped_reads_bwrd table.
    # 1: left matrix
    counts_l = _resore_orig_counts(counts_l, original_rc_l, pos_l)

    # 2: right matrix
    counts_r = _resore_orig_counts(counts_r, original_rc_r, pos_r)

    # Collect numbers of reads at positions for left junctions.
    pairs_df['N_unclipped_l'] = pairs_df.apply(
        lambda pos: unclipped_depth[pos.Chrom].get(pos.Position_l, 0),
        axis=1
    )
    pairs_df['N_clipped_l'] = pairs_df.apply(
        lambda pair: counts_l[
            pos_l[pair.Chrom][pair.Position_l], is_names_l[pair.IS_name]] if pair.Position_l > 0 else 0,
        axis=1
    )

    # Collect numbers of reads at positions for right junctions.
    pairs_df['N_unclipped_r'] = pairs_df.apply(
        lambda pos: unclipped_depth[pos.Chrom].get(pos.Position_r, 0),
        axis=1
    )
    pairs_df['N_clipped_r'] = pairs_df.apply(
        lambda pair: counts_r[
            pos_r[pair.Chrom][pair.Position_r], is_names_r[pair.IS_name]
        ] if pair.Position_r > 0 else 0,
        axis=1
    )

    # Add coverage from clipped reads that overlap with junction.
    pairs_df['N_overlap_l'] = pairs_df[['Position_l', 'Chrom']]. \
        apply(lambda x: cl_read_cov_overlap[x.Chrom].get(x.Position_l, 0), axis=1)
    pairs_df['N_overlap_r'] = pairs_df[['Position_r', 'Chrom']]. \
        apply(lambda x: cl_read_cov_overlap[x.Chrom].get(x.Position_r, 0), axis=1)

    # Metrics for corrections and tests
    min_match = min(match_lengths)
    av_read_len = read_lengths / n_reads_analyzed

    # Add corrections for clipped reads.
    pairs_df['N_clipped_l_correction'] = pairs_df['N_clipped_l'] / \
                                          (1 - min_match / av_read_len) - \
                                          pairs_df['N_clipped_l']

    pairs_df['N_clipped_r_correction'] = pairs_df['N_clipped_r'] / \
                                          (1 - min_match / av_read_len) - \
                                          pairs_df['N_clipped_r']

    pairs_df['N_overlap_l_correction'] = pairs_df['N_overlap_l'] / \
                                          (1 - min_match / av_read_len) - \
                                          pairs_df['N_overlap_l']

    pairs_df['N_overlap_r_correction'] = pairs_df['N_overlap_r'] / \
                                          (1 - min_match / av_read_len) - \
                                          pairs_df['N_overlap_r']

    pairs_df['N_clipped_l_corrected'] = pairs_df['N_clipped_l'] + pairs_df['N_clipped_l_correction']
    pairs_df['N_overlap_l_corrected'] = pairs_df['N_overlap_l'] + pairs_df['N_overlap_l_correction']
    pairs_df['N_clipped_r_corrected'] = pairs_df['N_clipped_r'] + pairs_df['N_clipped_r_correction']
    pairs_df['N_overlap_r_corrected'] = pairs_df['N_overlap_r'] + pairs_df['N_overlap_r_correction']

    pairs_df['N_overlap_formula_l'] = pairs_df[
        ['N_overlap_l_corrected', 'N_clipped_r_corrected', 'Position_l', 'Position_r']
    ]. \
        apply(lambda x: x.N_overlap_l_corrected - x.N_clipped_r_corrected
                        if x.Position_r > x.Position_l > 0
                        else x.N_overlap_l_corrected,
              axis=1)

    pairs_df['N_overlap_formula_r'] = pairs_df[
        ['N_overlap_r_corrected', 'N_clipped_l_corrected', 'Position_l', 'Position_r']
    ]. \
        apply(lambda x: x.N_overlap_r_corrected - x.N_clipped_l_corrected
                        if x.Position_r > x.Position_l > 0
                        else x.N_overlap_r_corrected,
              axis=1)

    # Calculate frequency as average between left and right boundaries if present.
    # If not - just by one boundary.
    # 0.1 pseudocount keeps from div/0 error.
    pairs_df['Frequency_l'] = pairs_df['N_clipped_l_corrected'] / \
                               (pairs_df['N_unclipped_l'] +
                                pairs_df['N_overlap_formula_l'] +
                                pairs_df['N_clipped_l_corrected'] +
                                0.1
                                )
    pairs_df['Frequency_r'] = pairs_df['N_clipped_r_corrected'] / \
                               (pairs_df['N_unclipped_r'] +
                                pairs_df['N_overlap_formula_r'] +
                                pairs_df['N_clipped_r_corrected'] +
                                0.1
                                )

    pairs_df['Frequency'] = pairs_df[['Frequency_l', 'Frequency_r']]. \
        apply(lambda x: _calc_freq_precise(x[0], x[1]), axis=1)

    # Add total depth column.
    pairs_df['Depth'] = pairs_df. \
        apply(
        lambda event:
        _add_total_depth(
            event.N_unclipped_l + event.N_overlap_formula_l + event.N_clipped_l_corrected,
            event.N_unclipped_r + event.N_overlap_formula_r + event.N_clipped_r_corrected
        ),
        axis=1
    )

    return pairs_df
