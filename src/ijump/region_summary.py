# Average mode's report generation: summarize junction counts by annotated
# region, then turn that summary into a per-region insertion-frequency report.

import logging

import pandas as pd


# Create an empty summary-by-region dataframe with one column per IS element.
def _sum_by_reg_tbl_init(is_coords):
    sbrcolumns = ['ann', 'chrom', 'start', 'stop']
    sbrcolumns.extend(list(is_coords.keys()))
    return pd.DataFrame(columns=sbrcolumns)


# Make summary table: for each annotated region, count supporting junctions
# per IS element.
def summarize_by_region(junctions, is_coords, gff_ann_pos) -> pd.DataFrame:
    logging.info('Create summary table by region')
    sum_by_region = pd.DataFrame()
    junc_temp = junctions.loc[junctions['Note'] != 'IS element']
    f_columns = ['ann', 'chrom', 'start', 'stop']
    f_columns.extend(list(is_coords.keys()))
    for i in range(len(junc_temp)):
        pos = junc_temp.iloc[i]['Position']
        chrom = junc_temp.iloc[i]['Chrom']
        for ann_id, item in gff_ann_pos[chrom].items():  #
            if item[2] <= pos <= item[3]:
                if ann_id not in sum_by_region.index:
                    columns = ['ann_id', 'ann', 'chrom', 'start', 'stop']
                    columns.extend(list(is_coords.keys()))
                    temp = _sum_by_reg_tbl_init(is_coords)
                    temp.at[0, 'ann_id'] = ann_id
                    temp.at[0, 'ann'] = item[0]
                    temp.at[0, 'chrom'] = item[1]
                    temp.at[0, 'start'] = item[2]
                    temp.at[0, 'stop'] = item[3]
                    temp.loc[0, list(is_coords.keys())] = 0
                    temp.at[0, junc_temp.iloc[i]['IS name']] = 1
                    temp = temp.set_index('ann_id')
                    sum_by_region = pd.concat([sum_by_region, temp], sort=True)
                else:
                    is_name = junc_temp.iloc[i]['IS name']
                    sum_by_region.loc[ann_id, is_name] += 1
                break
    sum_by_region = sum_by_region[f_columns]
    return sum_by_region


# Create report by IS and region: melt the region summary into one row per
# (IS, region) insertion event and add depth/frequency columns.
def report_average(sum_by_region, match_lengths, read_lengths, n_reads_analyzed,
                    blast_min, min_match, average_depth) -> pd.DataFrame:
    logging.info("Create report table")
    min_match = min(match_lengths)  # find minimum match length
    av_read_len = read_lengths / n_reads_analyzed  # find average read length
    report_table = pd.melt(
        sum_by_region,
        id_vars=('ann', 'chrom', 'start', 'stop'),
        var_name='IS Name',
        value_name='count'
    )

    # Drop zero intervals.
    report_table['drop'] = report_table.apply(lambda x: 0 if x['stop'] - x['start'] > 0 else 1, axis=1)
    report_table = report_table[report_table['drop'] == 0]
    report_table.drop(columns='drop', inplace=True)
    report_table.sort_values(by=['start', 'stop'], inplace=True)
    report_table = report_table[report_table['count'] > 0]

    # Add depth.
    report_table['Depth'] = report_table.apply(
        lambda x: average_depth(x['chrom'], x['start'], x['stop']),
        axis=1
    )

    report_table['Frequency'] = report_table.apply(
        lambda x: round((x['count'] / 2 * (1 + blast_min / av_read_len)) / (
                x['Depth'] * (1 - min_match / av_read_len)), 4),
        axis=1
    )

    report_table = report_table[['IS Name', 'ann', 'chrom', 'start', 'stop', 'Frequency', 'Depth']]
    report_table.columns = ['IS Name', 'Annotation', 'Chromosome', 'Start', 'Stop', 'Frequency', 'Depth']

    return report_table
