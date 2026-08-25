# Circos file generation: renders already-finished detection results (report_table,
# sum_by_region, is_coords, ref_len) into the file set an external tool (Circos) reads.
# Unrelated to IS-rearrangement detection itself -- see ticket 11.

import logging
import os
import re


# Create Circos files.
#
# `average_depth_fn` and `gff_ann_pos` are not part of ticket 11's originally enumerated
# signature (report_table, sum_by_region, is_coords, ref_len, data_folder, cutoff) but the
# moved body reads both: average-depth lookups (originally ISClipped.average_depth, which
# needs a live BAM handle via pysam) and GFF annotation positions (originally
# self.gff.ann_pos, a plain nested dict). Both are passed in explicitly instead of an object
# to read from, keeping this a plain function of its inputs -- true to the "already-finished
# results" framing in the ticket's Why section, just with two more of them than were
# enumerated there. `average_depth_fn` is a callable (chrom, start, stop) -> depth, not a
# scalar depth value -- named with the `_fn` suffix to keep that distinct from a plain
# average-depth number.
def write_files(
    report_table,
    sum_by_region,
    is_coords,
    ref_len,
    data_folder,
    cutoff,
    average_depth_fn,
    gff_ann_pos,
) -> None:
    logging.info("Create CIRCOS files")
    while not os.path.exists(data_folder):
        os.makedirs(data_folder)

    # Colors used for IS elements and contig representation
    _circos_colors = ("green", "red", "blue", "purple", "orange", "yellow", "grey")
    # Colors assigned to each chromosome
    _ref_colours = dict()
    # Colours assigned to each IS element
    _is_colours = dict()

    # Karyotype file
    with open(os.path.join(data_folder, "karyotype.txt"), "w") as karyotype:
        col_ind = 0
        for contig in ref_len.keys():
            karyotype.write(
                "chr - "
                + contig
                + " "
                + contig
                + " 0 "
                + str(ref_len[contig])
                + " "
                + _circos_colors[col_ind % len(_circos_colors)]
                + "\n"
            )
            _ref_colours[contig] = _circos_colors[col_ind % len(_circos_colors)]
            col_ind += 1

    # Text file
    with open(os.path.join(data_folder, "text.txt"), "w") as text:
        col_ind = 0
        for is_name in is_coords.keys():
            text.write(
                is_coords[is_name][0]
                + " "
                + is_coords[is_name][1]
                + " "
                + is_coords[is_name][1]
                + " "
                + is_name
                + " color=vvd"
                + _circos_colors[col_ind % len(_circos_colors)]
                + "\n"
            )
            _is_colours[is_name] = _circos_colors[col_ind % len(_circos_colors)]
            col_ind += 1

        # List to remove duplicates
        text_regions = list()

        # Add regions information.
        for i in range(len(report_table)):
            # Draw only lines with cutoff more than specified.
            if report_table.iloc[i]["Frequency"] >= cutoff:
                chrom = report_table.iloc[i]["Chromosome"]
                pos = report_table.iloc[i]["Start"]
                ann = report_table.iloc[i]["Annotation"]
                for a in ann.split("<>"):
                    if a in text_regions:
                        continue
                    text_regions.append(a)
                    text.write(chrom + " " + str(pos) + " " + str(pos) + " " + a + "\n")

    # links
    with open(data_folder + "links.txt", "w") as links:
        for i in range(len(report_table)):
            # Draw only lines with cutoff more than specified
            if report_table.iloc[i]["Frequency"] >= cutoff:
                is_name = report_table.iloc[i]["IS Name"]
                # Slice rather than unpack the whole entry: the IS table it
                # comes from grows columns over time, and Circos wants the
                # coordinate triple, not everything the table knows.
                is_chrom, is_start, is_stop = is_coords[is_name][:3]
                j_chrom = report_table.iloc[i]["Chromosome"]
                j_pos = str(report_table.iloc[i]["Start"])
                colour = "l" + _is_colours[is_name]
                links.write(
                    is_chrom
                    + " "
                    + is_start
                    + " "
                    + is_stop
                    + " "
                    + j_chrom
                    + " "
                    + j_pos
                    + " "
                    + j_pos
                    + " color="
                    + colour
                    + "\n"
                )

    # Histogram file
    with open(data_folder + "histogram.txt", "w") as histogram:
        for i in range(len(sum_by_region)):
            # Calculate average depth of the region.
            depth = average_depth_fn(
                sum_by_region.iloc[i]["chrom"],
                sum_by_region.iloc[i]["start"],
                sum_by_region.iloc[i]["stop"],
            )

            # Recalculate junction counts to depth.
            h_columns = ["chrom", "start", "stop"]
            h_columns_is = [x for x in is_coords.keys()]

            for h in h_columns_is:
                if depth > 0:
                    if sum_by_region.iloc[i][h] / depth / 2 >= cutoff:
                        histogram.write(
                            " ".join(sum_by_region.iloc[i][h_columns].apply(str))
                            + " "
                            + ",".join(
                                sum_by_region.iloc[i][h_columns_is]
                                .apply(lambda x, depth=depth: round(((x / depth / 2) * 100), 2))
                                .apply(str)
                            )
                            + "\n"
                        )
                        break

    # Depth histogram
    with open(data_folder + "depth.txt", "w") as depth_hist:
        for contig in gff_ann_pos:
            for _ann_id, ann in gff_ann_pos[contig].items():
                if ann[3] - ann[2] <= 0:
                    continue
                depth = average_depth_fn(ann[1], ann[2], ann[3])
                depth_hist.write(" ".join([str(x) for x in ann[1:]]) + " " + str(depth) + "\n")

    # Write config.
    config_name = os.path.join(data_folder, "circos.conf")
    with open(config_name, "w") as config:
        # circos.conf lives at the repo root; this module lives two levels
        # below it at src/ijump/circos.py.
        repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
        logging.info(repo_root)
        conf_template = open(repo_root + "/circos.conf", "r")
        conf = conf_template.read()
        conf = re.sub("karyotype = XXX", "karyotype = " + data_folder + "karyotype.txt", conf)
        conf = re.sub("XXX		#text", data_folder + "text.txt", conf)
        conf = re.sub("XXX		#links", data_folder + "links.txt", conf)
        conf = re.sub("XXX		#histogram", data_folder + "histogram.txt", conf)
        conf = re.sub("XXX		#depth", data_folder + "depth.txt", conf)

        # Make fill_color string for a histogram.
        hist_colors = ""
        for is_name in is_coords.keys():
            hist_colors += _is_colours[is_name] + ", "
        hist_colors = hist_colors[:-2]

        conf = re.sub("XXX		#stacked_colors", hist_colors, conf)

        config.write(conf)
