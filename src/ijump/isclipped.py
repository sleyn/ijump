# Implement resolve position conflict function when two IS elements are inserted at
# the same position.

import logging
import os
import re
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from statistics import mean

import numpy as np
import pandas as pd
import pysamstats
from sklearn.cluster import AgglomerativeClustering

from ijump import circos, clipped_read_search, frequency_estimation, gff, region_summary
from ijump.clipped_read_search import NoInsertionsFound
from ijump.junction_pairing import find_pairs


# Values accepted by --estimation_mode. A str Enum instead of bare string
# literals so a mistyped comparison raises AttributeError/NameError instead
# of silently falling through a dead branch (see ijump_junctions.txt 'IS pos'
# regression fixed alongside this).
class EstimationMode(str, Enum):
    AVERAGE = "average"
    PRECISE = "precise"


# Result of ISClipped.run(): whether the pipeline found anything to report,
# and (if not) the reason, so the caller can log it without inspecting
# internals.
@dataclass
class RunResult:
    insertions_found: bool
    message: str = ""


# Check if any clipped reads are present.
def check_data_presence_in_df(cl_table, message):
    if cl_table.size == 0:
        raise NoInsertionsFound(message)


# Check if junctions exist for IS elements.
def check_junctions_presence(junc_tbl, outdir, est_mode):
    if junc_tbl.size:
        # Convert from 0-base to 1-base system
        junc_tbl_copy = junc_tbl.copy()
        if est_mode == EstimationMode.PRECISE:
            junc_tbl_copy["IS pos"] += 1
        junc_tbl_copy["Position"] += 1
        junc_tbl_copy.to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep="\t", index=False)
    else:
        raise NoInsertionsFound("No junctions was found")


# Claculate distance between insertion positions.
def interpos_distance(pos_l, pos_r):
    if pos_l == 0:
        return 5
    elif pos_r == 0:
        return 5
    else:
        return abs(pos_r - pos_l) + 5


# Assign Keep status to one pair.
def keep_pair(pair, region_starts, region_ends):
    for position in pair:
        compare_starts = region_starts <= position
        compare_ends = region_ends >= position
        if np.any(np.all([compare_starts, compare_ends], axis=0)):
            return True
    return False


# Filter pairs. Keep a pair if at least one position in pair is in the region interval.
def filter_pairs(pairs_tbl, region_tbl):
    logging.info("Filter pairs that do not belong to the regions of interest.")
    regions_starts = region_tbl["Position_left"].to_numpy()
    regions_ends = region_tbl["Position_right"].to_numpy()
    pairs_tbl_return = pairs_tbl.copy()
    pairs_tbl_return["Keep"] = pairs_tbl_return[["Position_l", "Position_r"]].apply(
        lambda pos_pair: keep_pair(pos_pair.to_list(), regions_starts, regions_ends), axis=1
    )
    return pairs_tbl_return[pairs_tbl_return["Keep"]].drop(columns=["Keep"])


# Convert coordinate system of a list from 0-base to 1-base.
def convert_zero_one_base(coordinates_column):
    return list(map(lambda x: x + 1 if x > 0 else 0, coordinates_column))


# specify class for clipped reads
class ISClipped:
    def __init__(
        self, aln, ref_name, gff_name, workdir="ijump_wd", pairs_df_path="ijump_junction_pairs.txt"
    ):
        # Input files:
        # pysam file
        self.aln = aln
        # File with the reference genome in FASTA format
        self.ref_name = ref_name
        # A file with GFF annotations
        self.gff_name = gff_name
        # Set a GFF object for GFF file.
        self.gff = gff.gff(self.gff_name)
        # Work directory directory
        self.workdir = workdir
        # Path to the output pairs table.
        # In case of preliminary exit an empty table will be generated.
        self.pairs_df_path = pairs_df_path

        # Tables:
        # Clipped reads from the forward search (IS -> Ref_genome)
        self.clipped_reads = self._cltbl_init()
        # Clipped reads from the backward search for precise pipeline (Ref_genome -> IS)
        self.clipped_reads_bwrd = self._cltbl_init()
        # Filtered BLAST output from search juction location for unaligned part of clipped reads
        self.blastout_filtered = self._blastout_filtered_init()
        # Junction positions table
        self.junctions = self._jtbl_init()
        # Table with Region-centric representation: for each Region number of
        # supporting reads for each IS is shown
        self.sum_by_region = pd.DataFrame()
        # Table of each insertion event
        self.report_table = pd.DataFrame()
        # Create pairs table.
        # Used in a precise mode to collect information of insertion coordinates in reference.
        self.pairs_df = self._pairs_table_init()

        # Helper structures:
        # Length of reference contigs
        self.ref_len = dict()
        # Dictionary of depth attriduted to unclipped reads in given positions
        # unclipped reads do not contain "S" in CIGAR string
        self.unclipped_depth = {}
        # Depth of coverage from clipped reads on the position of junctions.
        # Only clipped reads that DO NOT have junction in the position of interest are count.
        # Required to separate close insertions of IS elements.
        # Populated by clipped_read_search.search's backward (direction=0)
        # call -- see SearchResult.cl_read_cov_overlap.
        self.cl_read_cov_overlap = {}
        # Boundaries of clipped reads
        self.boundaries = list()
        # Coordinates of provided IS elements. IS name => [chrom, start, stop]
        self.is_coords = dict()
        # List of lengths for matched segments
        # Used to calculate correction coefficients
        # Populated by clipped_read_search.search's forward (direction=1) call.
        self.match_lengths = list()

        # Initialize nested dictionaries as for each dictionary we need contig->position id.
        for contig_i in range(len(self.aln.references)):
            self.ref_len[self.aln.references[contig_i]] = self.aln.lengths[contig_i]
            self.unclipped_depth[self.aln.references[contig_i]] = {}

        # Parameters:
        # Show junctions only with this frequency or more
        self.cutoff = 0.005
        # Average read length
        self.av_read_len = 150
        # Total length of reads
        # Populated by clipped_read_search.search's forward (direction=1) call.
        self.read_lengths = 0
        # Number of analyzed reads
        # Populated by clipped_read_search.search's forward (direction=1) call.
        self.n_reads_analyzed = 0
        # Minimum length of clipped part to use in BLAST. Fed into the
        # frequency calculation in region_summary.report_average; the value
        # itself lives on clipped_read_search.BLAST_MIN.
        self.blast_min = clipped_read_search.BLAST_MIN
        # Maximum expected length of duplication created from the insertion event
        self.max_is_dup_len = 20

        # Data folder for circos files
        self.data_folder = "./ijump_data/"

    # Initialize a pairs table.
    # Used in a precise mode to collect information of insertion coordinates in reference.
    @staticmethod
    def _pairs_table_init():
        return pd.DataFrame(  # prototype of pairs table
            {
                "IS_name": ["-"],
                "Position_l": [0],
                "Position_r": [0],
                "Count_mapped_to_IS_l": [0],
                "Count_mapped_to_IS_r": [0],
                "Chrom": ["-"],
            }
        )

    # Generate an empty output table
    @staticmethod
    def pairs_table_empty():
        return pd.DataFrame(
            {
                "Position_l": [],
                "Position_r": [],
                "Count_mapped_to_IS_l": [],
                "Count_mapped_to_IS_r": [],
                "Chrom": [],
                "IS_name": [],
                "Dist": [],
                "N_unclipped_l": [],
                "N_clipped_l": [],
                "N_unclipped_r": [],
                "N_clipped_r": [],
                "N_overlap_l": [],
                "N_overlap_r": [],
                "N_clipped_l_correction": [],
                "N_clipped_r_correction": [],
                "N_overlap_l_correction": [],
                "N_overlap_r_correction": [],
                "N_clipped_l_corrected": [],
                "N_overlap_l_corrected": [],
                "N_clipped_r_corrected": [],
                "N_overlap_r_corrected": [],
                "N_overlap_formula_l": [],
                "N_overlap_formula_r": [],
                "Frequency_l": [],
                "Frequency_r": [],
                "Frequency": [],
                "Depth": [],
            }
        )

    # Initialize a report table.
    @staticmethod
    def report_table_init():
        return pd.DataFrame(
            columns=["IS Name", "Annotation", "Chromosome", "Start", "Stop", "Frequency", "Depth"]
        )

    # Initialize a new copy of clipped reads table.
    @staticmethod
    def _cltbl_init():
        return pd.DataFrame(
            columns=[
                "ID",
                "IS name",
                "IS_chrom",
                "Read name",
                "left pos",  # left position of a clipped segment in a read
                "right pos",  # right position of a clipped segment in a read
                "clip_position",  # clip position in a reference
                "junction_in_read",  # side of a clipped segment connected to a read (l/r)
                "reverse",  # is read reverse
                "sequence",
            ]
        )  # sequence of a clipped read

    # Initialize a filtered blast output table.
    @staticmethod
    def _blastout_filtered_init():
        return pd.DataFrame(
            columns=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "pos_in_ref",
                "orientation",
            ]
        )

    # Create summary dataframe. Shape is shared with region_summary's
    # per-region builder -- see region_summary.sum_by_reg_tbl_init.
    def sum_by_reg_tbl_init(self):
        return region_summary.sum_by_reg_tbl_init(self.is_coords)

    # Collect information about IS elements.
    def iscollect(self, file):
        logging.info(f"Read file with IS elements: {file}")
        with open(file, "r") as is_coords_file:
            for coord in is_coords_file.readlines():
                c = coord.split()
                self.is_coords[c[0]] = c[1:]

    # Initialize junction table.
    @staticmethod
    def _jtbl_init(n_rows=0):
        return pd.DataFrame(
            columns=[
                "ID",
                "IS name",
                "IS pos",
                "IS_chrom",
                "Read name",
                "Chrom",
                "Position",
                "Orientation",  # where non-clipped region is relative to position
                "Note",
                "Locus tag",
                "Gene",
            ],
            index=[i for i in range(n_rows)],
        )

    # For each IS element set boundaries where to search clipped reads.
    def set_is_boundaries(self, radius):
        logging.info("Set area near IS elements boundaries to search clipped reads.")
        for is_name in self.is_coords.keys():
            # For a IS element take its contig and coordinates.
            chrom = self.is_coords[is_name][0]
            start = int(self.is_coords[is_name][1])
            stop = int(self.is_coords[is_name][2])

            # Check if search area for clipped reads goes outside contig boundaries.
            # ASSUMPTION OF COMPLETENESS!: If it does so consider all contigs as circular.
            # Check for a start coordinate.
            if start - radius < 0:
                self.boundaries.append(
                    [
                        self.ref_len[chrom] - radius + start,
                        self.ref_len[chrom],
                        "start",
                        is_name,
                        chrom,
                    ]
                )
                self.boundaries.append([0, start + radius, "start", is_name, chrom])
            elif start + radius > self.ref_len[chrom]:
                self.boundaries.append(
                    [start - radius, self.ref_len[chrom], "start", is_name, chrom]
                )
                self.boundaries.append(
                    [0, start + radius - self.ref_len[chrom], "start", is_name, chrom]
                )
            else:
                self.boundaries.append([start - radius, start + radius, "start", is_name, chrom])

            # Check for an end coordinate.
            if stop + radius > self.ref_len[chrom]:
                self.boundaries.append([stop - radius, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append(
                    [0, stop + radius - self.ref_len[chrom], "stop", is_name, chrom]
                )
            elif stop - radius < 0:
                self.boundaries.append(
                    [
                        self.ref_len[chrom] - radius + stop,
                        self.ref_len[chrom],
                        "stop",
                        is_name,
                        chrom,
                    ]
                )
                self.boundaries.append([0, stop + radius, "stop", is_name, chrom])
            else:
                self.boundaries.append([stop - radius, stop + radius, "stop", is_name, chrom])

    # Check if position close to the IS element
    def _check_is_boundary_proximity(self, chrom, position):
        for b in self.boundaries:
            if b[4] == chrom:
                # if b[0] - boundary_width / 2 <= position <= b[1] + boundary_width / 2:
                # use doubled boundaries
                if b[0] <= position <= b[1]:
                    return "IS element", b[3]
        return "-", "-"

    # Cluster positions together to form seed for backwards clipped read search.
    @staticmethod
    def _hclust(X):
        # If only one sample is present – clustering will not work
        if len(X) == 1:
            return [0]
        hcl = AgglomerativeClustering(n_clusters=None, distance_threshold=30, linkage="single").fit(
            X.to_numpy().reshape(-1, 1)
        )
        return hcl.labels_

    # Use hierarchical clustering to cluster close positions in the chromosome.
    def make_gene_side_regions(self):
        logging.info("Cluster close positions in the chromosome")

        # Remove positions close to the IS elements boundaries from the analysis
        ref_cl_reads = self.blastout_filtered[["sseqid", "pos_in_ref"]].copy()
        ref_cl_reads = ref_cl_reads.rename(columns={"sseqid": "Chrom", "pos_in_ref": "Position"})
        ref_cl_reads["Note"] = ref_cl_reads.apply(
            lambda x: self._check_is_boundary_proximity(x["Chrom"], x["Position"])[0], axis=1
        )
        ref_cl_reads = ref_cl_reads[ref_cl_reads["Note"] == "-"]
        ref_cl_reads = ref_cl_reads.drop(columns=["Note"])
        # If no hits point outside IS elements boudaries there is no insertions to find
        if ref_cl_reads.size == 0:
            logging.info(
                "No BLAST hits point oustide IS elements. No significant new insertions"
                " could be found."
            )
            raise NoInsertionsFound(
                "No BLAST hits point oustide IS elements. No significant new insertions"
                " could be found."
            )

        ref_cl_reads["Cluster"] = (
            ref_cl_reads.sort_values(by=["Chrom", "Position"])
            .groupby(["Chrom"])["Position"]
            .transform(self._hclust)
        )

        ref_regions = (
            ref_cl_reads.groupby(["Cluster", "Chrom"]).aggregate(["min", "max"]).reset_index()
        )

        ref_regions.columns = ["Cluster", "Chrom", "Position_left", "Position_right"]
        ref_regions = ref_regions.drop(columns=["Cluster"])

        # Extend regions by 5nt if possible
        ref_regions["Position_left"] = ref_regions["Position_left"].apply(lambda x: max(x - 5, 0))

        ref_regions["Position_right"] = ref_regions.apply(
            lambda x: min(x["Position_right"] + 5, self.ref_len[x["Chrom"]]), axis=1
        )

        return ref_regions

    # Create table for description of junctions.
    # direction: 1=>IS->Ref, 0=>Ref->IS.
    def call_junctions(self, direction):
        logging.info("Create junction table")
        self.junctions = self._jtbl_init(self.blastout_filtered.shape[0])
        index = 0

        for hit in self.blastout_filtered.itertuples(index=False):
            read_id = hit.qseqid
            if direction:
                pos = hit.pos_in_ref
                chrom = hit.sseqid
                is_name = self.clipped_reads["IS name"][read_id]
                is_chrom = self.clipped_reads["IS_chrom"][read_id]
                is_pos = self.clipped_reads["clip_position"][read_id]
                is_elem_border_mark, _ = self._check_is_boundary_proximity(chrom, pos)
                orientation = hit.orientation
                read_name = self.clipped_reads["Read name"][read_id]
            else:
                pos = self.clipped_reads_bwrd["junction_in_read"][read_id]
                chrom = self.clipped_reads_bwrd["IS_chrom"][read_id]
                is_chrom = hit.sseqid
                is_pos = hit.pos_in_ref
                _, is_name = self._check_is_boundary_proximity(is_chrom, is_pos)
                is_elem_border_mark = "-"
                orientation = self.clipped_reads_bwrd["clip_position"][read_id]
                read_name = self.clipped_reads_bwrd["Read name"][read_id]

            self.junctions.at[index, "ID"] = read_id
            self.junctions.at[index, "IS name"] = is_name
            self.junctions.at[index, "IS_chrom"] = is_chrom
            self.junctions.at[index, "IS pos"] = is_pos
            self.junctions.at[index, "Read name"] = read_name
            self.junctions.at[index, "Chrom"] = chrom
            self.junctions.at[index, "Position"] = pos
            self.junctions.at[index, "Orientation"] = orientation
            self.junctions.at[index, "Locus tag"] = self.gff.gff_pos[chrom][pos][0]
            self.junctions.at[index, "Gene"] = self.gff.gff_pos[chrom][pos][1]
            self.junctions.at[index, "Note"] = is_elem_border_mark
            index += 1

        self.junctions = self.junctions.reset_index()

    # Find positions of insertions.
    def search_insert_pos(self):
        logging.info("Serach for junction pairs")
        position_tbl = self.junctions[self.junctions["IS name"] != "-"].copy()
        # It is much better to work with when IS elements collapsed by their copy
        # than to work with each copy separately.
        # Remove copy tags from the IS element names like "_1", "_2".
        position_tbl["IS"] = position_tbl["IS name"].apply(
            lambda x: re.search(r"(.+)_\d+", x).group(1)
        )
        position_tbl = (
            position_tbl.groupby(["Chrom", "Position", "IS", "Orientation"])["Position"]
            .count()
            .reset_index(name="Counts")
        )

        # Collect dataframes for pairs of junctions (or orphan junctions) that should
        # mark IS elements insertions.
        is_pairs_collection = []

        # Find pairs
        for chrom in position_tbl["Chrom"].drop_duplicates().tolist():
            # Take IS elements only from the selected chromosome.
            position_tbl_chrom = position_tbl.query("Chrom == @chrom")
            for is_name in position_tbl_chrom["IS"].drop_duplicates().tolist():
                positions_left = position_tbl_chrom.query(
                    'Orientation == "left" & IS == @is_name'
                ).sort_values("Position")
                positions_left_pos = positions_left["Position"].to_numpy()
                positions_left_counts = positions_left["Counts"].to_numpy()
                positions_right = position_tbl_chrom.query(
                    'Orientation == "right" & IS == @is_name'
                ).sort_values("Position")
                positions_right_pos = positions_right["Position"].to_numpy()
                positions_right_counts = positions_right["Counts"].to_numpy()

                logging.info(f"Find pairs for {is_name} and {chrom} contig ")
                # Calculate table of pairs
                pair_tbl_chunk = find_pairs(
                    positions_left_pos,
                    positions_right_pos,
                    positions_left_counts,
                    positions_right_counts,
                    self.ref_len[chrom],
                    self.max_is_dup_len,
                    chrom,
                )

                pair_tbl_chunk["IS_name"] = is_name

                is_pairs_collection.append(pair_tbl_chunk)

        # Concatenate all pair tables into one table.
        self.pairs_df = pd.concat(is_pairs_collection, ignore_index=True)

    # Count depth at the position using only unclipped reads.
    # Input is a data frame with columns 'Position' and 'Chrom'.
    def count_depth_unclipped(self, position_tbl):
        logging.info("Count depth attributed to unclipped reads")
        for position in position_tbl.itertuples():
            chrom = position.Chrom
            pos = position.Position
            ins_pos_distance = position.Dist

            if pos == 0:
                continue

            for read in self.aln.fetch(chrom, pos, pos + 1):
                if read.is_unmapped:
                    # Skip unmapped read
                    continue
                # No soft- or hard-clipped reads
                elif ("S" not in read.cigarstring) and ("H" not in read.cigarstring):
                    # Test if position is near the end of the read. If it is near skip the read as
                    # it is impossible to distinguish unmapped
                    read_edges = np.array([read.aligned_pairs[0][1], read.aligned_pairs[-1][1]])
                    if np.min(np.abs(pos - read_edges)) > ins_pos_distance:
                        self.unclipped_depth[chrom][pos] = (
                            self.unclipped_depth[chrom].get(pos, 0) + 1
                        )

    # Write the full set of output files with headers and zero data rows.
    # A run that finds nothing is a successful run, so it writes the same
    # file set a run that finds something would.
    def _write_empty_outputs(self, mode, outdir):
        self._cltbl_init().to_csv(os.path.join(outdir, "reads.txt"), sep="\t", index=False)
        self._jtbl_init().to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep="\t", index=False)
        self.sum_by_reg_tbl_init().to_csv(
            os.path.join(outdir, "ijump_sum_by_reg.txt"), sep="\t", index=False
        )
        self.report_table_init().to_csv(
            os.path.join(outdir, "ijump_report_by_is_reg.txt"), sep="\t", index=False
        )

        if mode == EstimationMode.PRECISE:
            self.pairs_table_empty().to_csv(
                os.path.join(outdir, "ijump_junction_pairs.txt"), sep="\t", index=False
            )

    # Run the full average/precise pipeline and write every output file it
    # produces along the way. A run that finds nothing to report is not an
    # error: it is signalled through the returned RunResult rather than by
    # letting NoInsertionsFound cross this method's boundary.
    def run(self, mode):
        outdir = os.path.dirname(self.pairs_df_path)

        try:
            # Collect clipped reads and search insertion positions in Reference.
            # 1 - search in IS->Reference direction.
            result = clipped_read_search.search(
                self.aln, self.boundaries, self.ref_name, self.workdir, direction=1
            )
            self.clipped_reads = result.clipped_reads
            self.blastout_filtered = result.blast_hits
            self.match_lengths = result.match_lengths
            self.read_lengths = result.read_lengths
            self.n_reads_analyzed = result.n_reads_analyzed
            self.clipped_reads.to_csv(os.path.join(outdir, "reads.txt"), sep="\t", index=False)

            # Read GFF file.
            self.gff.readgff()

            if mode == EstimationMode.AVERAGE:
                # Make a table of observed junction positions
                self.call_junctions(1)

                # Check if any junction is present. If not - stop the workflow.
                check_junctions_presence(self.junctions, outdir, mode)

                # Count reads supporting IS elements insertions for each IS element and each Region
                # Reformat GFF representation
                self.gff.pos_to_ann()
                self.sum_by_region = region_summary.summarize_by_region(
                    self.junctions, self.is_coords, self.gff.ann_pos
                )
                self.sum_by_region.to_csv(
                    os.path.join(outdir, "ijump_sum_by_reg.txt"), sep="\t", index=False
                )

                # Make a report table of assessed insertion frequencies in each Region
                self.report_table = region_summary.report_average(
                    self.sum_by_region,
                    self.match_lengths,
                    self.read_lengths,
                    self.n_reads_analyzed,
                    self.blast_min,
                    self.average_depth,
                )
                self.report_table.to_csv(
                    os.path.join(outdir, "ijump_report_by_is_reg.txt"), sep="\t", index=False
                )
            elif mode == EstimationMode.PRECISE:
                # Make table of regions in the reference genome where extract clipped
                # reads for backwards assignment
                reference_regions = self.make_gene_side_regions()
                reference_regions.to_csv(
                    os.path.join(self.workdir, "reference_regions.tsv"), sep="\t", index=False
                )

                # Collect clipped reads at the insertion positions found during
                # forward (IS->Reference) search, then search for their
                # positions in the reference (Reference->IS direction).
                backward_boundaries = [
                    [row.Position_left, row.Position_right, "-", "-", row.Chrom]
                    for row in reference_regions.itertuples()
                ]
                result = clipped_read_search.search(
                    self.aln,
                    backward_boundaries,
                    self.ref_name,
                    self.workdir,
                    direction=0,
                )
                self.clipped_reads_bwrd = result.clipped_reads
                self.blastout_filtered = result.blast_hits
                self.cl_read_cov_overlap = result.cl_read_cov_overlap

                # Format results as junction table to attribute reads and their
                # junction positions to the IS elements.
                self.call_junctions(0)

                # Check if any junction is present. If not - stop the workflow.
                check_junctions_presence(self.junctions, outdir, mode)

                # Find pairs of junctions that should indicate insertion positions of
                # both edges of IS element.
                self.search_insert_pos()

                # Filter Junction pairs so at least one of the pair is in the
                # "reference_regions" table.
                self.pairs_df = filter_pairs(self.pairs_df, reference_regions)

                # Check if any pair was produced
                check_data_presence_in_df(self.pairs_df, "No pairs were found.")

                # Count depth of unclipped reads to have a background depth of coverage
                # Preparation.
                self.pairs_df["Dist"] = self.pairs_df[["Position_l", "Position_r"]].apply(
                    lambda clustered_pos: interpos_distance(
                        clustered_pos.Position_l, clustered_pos.Position_r
                    ),
                    axis=1,
                )

                positions = pd.concat(
                    [
                        self.pairs_df[["Position_l", "Chrom", "Dist"]].rename(
                            columns={"Position_l": "Position"}
                        ),
                        self.pairs_df[["Position_r", "Chrom", "Dist"]].rename(
                            columns={"Position_r": "Position"}
                        ),
                    ],
                    axis=0,
                ).drop_duplicates()

                # The depth count itself
                self.count_depth_unclipped(positions)
                self.pairs_df.drop(columns="Dist")

                # Make an estimate of insertion frequency
                logging.info("Estimate insertion frequencies.")
                self.pairs_df = frequency_estimation.estimate_frequencies(
                    self.pairs_df,
                    self.clipped_reads_bwrd,
                    self.unclipped_depth,
                    self.cl_read_cov_overlap,
                    self.match_lengths,
                    self.read_lengths,
                    self.n_reads_analyzed,
                )

                # Convert coordinates from 0-base to 1-base
                self.pairs_df["Position_l"] = convert_zero_one_base(
                    self.pairs_df["Position_l"].tolist()
                )
                self.pairs_df["Position_r"] = convert_zero_one_base(
                    self.pairs_df["Position_r"].tolist()
                )
                self.pairs_df.to_csv(self.pairs_df_path, sep="\t", index=False)
        except NoInsertionsFound as e:
            message = str(e)
            logging.info(message)
            self._write_empty_outputs(mode, outdir)
            return RunResult(insertions_found=False, message=message)

        return RunResult(insertions_found=True)

    # Calculate average depth of the region.
    # NOTE (ticket 06): caching on an instance method via `lru_cache` keeps a
    # reference to `self` for the cache's lifetime, which can leak memory if
    # many ISClipped instances are created. Flagged rather than restructured,
    # per ticket 06's scope (tooling-config only, no design changes).
    @lru_cache(maxsize=128)  # noqa: B019
    def average_depth(self, chrom, start, stop):
        # aln_depth = self.aln.count_coverage(chrom, start, stop)
        # depth = sum(map(sum, aln_depth))
        # return depth / len(aln_depth[0])  # average depth of the region
        c = pysamstats.load_coverage(
            self.aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000
        )
        return mean(c.reads_all)

    # Write the Circos diagram files from this pipeline's own state.
    # The caller (main()) only decides whether to call this; this method
    # is the adapter that knows which pieces of pipeline state a Circos
    # diagram needs.
    def write_circos_files(self):
        circos.write_files(
            self.report_table,
            self.sum_by_region,
            self.is_coords,
            self.ref_len,
            self.data_folder,
            self.cutoff,
            self.average_depth,
            self.gff.ann_pos,
        )
