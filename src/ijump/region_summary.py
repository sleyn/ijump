# Average mode's report generation: summarize junction counts by annotated
# region, then turn that summary into a per-region insertion-frequency report.

import logging

import numpy as np
import pandas as pd


# The columns a region summary carries, in order: the fixed region columns, then
# one per cluster.
#
# One column per *cluster*, not per locus. A read clipped at one locus of an
# element is indistinguishable from one clipped at another, so a column per locus
# splits one insertion's evidence across as many columns as the assembly happened
# to call fragments of that element -- three, on the reference genome, where
# IS17_1, IS17_2 and ISAba12_1 are one copy and two of its own fragments
# (isfinder-annotation 07).
#
# `clusters` maps each IS name to its cluster (is_table.cluster_by_name). Column
# order follows the IS table's row order, deduplicated.
def region_columns(clusters):
    return ["ann", "chrom", "start", "stop"] + list(dict.fromkeys(clusters.values()))


# Create an empty summary-by-region dataframe with one column per cluster.
# Shared by the no-results path (ISClipped._write_empty_outputs, via
# ISClipped.sum_by_reg_tbl_init) and the populated path (summarize_by_region
# below), so the two output files can never disagree on shape.
def sum_by_reg_tbl_init(clusters):
    return pd.DataFrame(columns=region_columns(clusters))


# Make summary table: for each annotated region, count supporting junctions
# per cluster.
def summarize_by_region(junctions, clusters, gff_ann_pos) -> pd.DataFrame:
    logging.info("Create summary table by region")
    sum_by_region = pd.DataFrame()
    junc_temp = junctions.loc[junctions["Note"] != "IS element"]
    f_columns = region_columns(clusters)
    cluster_names = f_columns[4:]
    for i in range(len(junc_temp)):
        pos = junc_temp.iloc[i]["Position"]
        chrom = junc_temp.iloc[i]["Chrom"]
        cluster = clusters[junc_temp.iloc[i]["IS name"]]
        for ann_id, item in gff_ann_pos[chrom].items():  #
            if item[2] <= pos <= item[3]:
                if ann_id not in sum_by_region.index:
                    temp = sum_by_reg_tbl_init(clusters)
                    temp.at[0, "ann_id"] = ann_id
                    temp.at[0, "ann"] = item[0]
                    temp.at[0, "chrom"] = item[1]
                    temp.at[0, "start"] = item[2]
                    temp.at[0, "stop"] = item[3]
                    temp.loc[0, cluster_names] = 0
                    temp.at[0, cluster] = 1
                    temp = temp.set_index("ann_id")
                    sum_by_region = pd.concat([sum_by_region, temp], sort=True)
                else:
                    sum_by_region.loc[ann_id, cluster] += 1
                break
    sum_by_region = sum_by_region[f_columns]
    return sum_by_region


# Create report by IS and region: melt the region summary into one row per
# (IS, region) insertion event and add depth/frequency columns.
def report_average(
    sum_by_region,
    match_lengths,
    read_lengths,
    n_reads_analyzed,
    blast_min,
    average_depth,
) -> pd.DataFrame:
    logging.info("Create report table")
    # 1st percentile rather than the bare minimum: match_lengths is a per-read
    # statistic pooled across the whole run, and a single outlier read (a
    # spuriously short matched segment) would otherwise set the correction
    # applied to every region's frequency (average-depth-zero-coverage 02).
    min_match = np.percentile(match_lengths, 1)
    av_read_len = read_lengths / n_reads_analyzed  # find average read length
    # "IS Name" carries the cluster -- the element that jumped. The column keeps
    # its name because it is the join key combine_results merges samples on and
    # the header of a documented output; what changed is that one name no longer
    # means one called locus.
    report_table = pd.melt(
        sum_by_region,
        id_vars=("ann", "chrom", "start", "stop"),
        var_name="IS Name",
        value_name="count",
    )

    # Drop zero intervals.
    report_table["drop"] = report_table.apply(
        lambda x: 0 if x["stop"] - x["start"] > 0 else 1, axis=1
    )
    report_table = report_table[report_table["drop"] == 0]
    report_table.drop(columns="drop", inplace=True)
    report_table.sort_values(by=["start", "stop"], inplace=True)
    report_table = report_table[report_table["count"] > 0]

    # Add depth.
    report_table["Depth"] = report_table.apply(
        lambda x: average_depth(x["chrom"], x["start"], x["stop"]), axis=1
    )

    report_table["Frequency"] = report_table.apply(
        lambda x: (
            np.nan
            if x["Depth"] == 0
            else round(
                (x["count"] / 2 * (1 + blast_min / av_read_len))
                / (x["Depth"] * (1 - min_match / av_read_len)),
                4,
            )
        ),
        axis=1,
    )

    report_table = report_table[["IS Name", "ann", "chrom", "start", "stop", "Frequency", "Depth"]]
    report_table.columns = [
        "IS Name",
        "Annotation",
        "Chromosome",
        "Start",
        "Stop",
        "Frequency",
        "Depth",
    ]

    return report_table
