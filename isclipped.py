# Implement resolve position conflict function when two IS elements are inserted at the same position.

import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
import re
import os
import gff
import frequency_estimation
from junction_pairing import find_pairs
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
import pysamstats
from statistics import mean
from sklearn.cluster import AgglomerativeClustering
import logging


class NoInsertionsFound(Exception):
    """The analysis ran correctly and there is nothing to report."""


# Values accepted by --estimation_mode. A str Enum instead of bare string
# literals so a mistyped comparison raises AttributeError/NameError instead
# of silently falling through a dead branch (see ijump_junctions.txt 'IS pos'
# regression fixed alongside this).
class EstimationMode(str, Enum):
    AVERAGE = 'average'
    PRECISE = 'precise'


# Result of ISClipped.run(): whether the pipeline found anything to report,
# and (if not) the reason, so the caller can log it without inspecting
# internals.
@dataclass
class RunResult:
    insertions_found: bool
    message: str = ''


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
            junc_tbl_copy['IS pos'] += 1
        junc_tbl_copy['Position'] += 1
        junc_tbl_copy.to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep='\t', index=False)
    else:
        raise NoInsertionsFound('No junctions was found')


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
    logging.info('Filter pairs that do not belong to the regions of interest.')
    regions_starts = region_tbl['Position_left'].to_numpy()
    regions_ends = region_tbl['Position_right'].to_numpy()
    pairs_tbl_return = pairs_tbl.copy()
    pairs_tbl_return['Keep'] = pairs_tbl_return[['Position_l', 'Position_r']].apply(
        lambda pos_pair: keep_pair(pos_pair.to_list(), regions_starts, regions_ends),
        axis=1
    )
    return pairs_tbl_return[pairs_tbl_return['Keep']].drop(columns=['Keep'])


# Convert coordinate system of a list from 0-base to 1-base.
def convert_zero_one_base(coordinates_column):
    return list(map(lambda x: x + 1 if x > 0 else 0, coordinates_column))


# specify class for clipped reads
class ISClipped:
    def __init__(self, aln, ref_name, gff_name, workdir='ijump_wd', pairs_df_path='ijump_junction_pairs.txt'):
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
        # Dictionary is used for speed purposes.
        self.clipped_reads = self._cltbl_init()
        self.clipped_reads_dict = {}
        # Clipped reads from the backward search for precise pipeline (Ref_genome -> IS)
        # Dictionary is used for speed purposes.
        self.clipped_reads_bwrd = self._cltbl_init()
        self.clipped_reads_bwrd_dict = {}
        # Filtered BLAST output from search juction location for unaligned part of clipped reads
        self.blastout_filtered = self._blastout_filtered_init()
        # Junction positions table
        self.junctions = self._jtbl_init()
        # Table with genetic element (GE) centic representation: for each GE number of supporting reads
        # for each IS is shown
        self.sum_by_region = pd.DataFrame()
        # Table of each insertion event
        self.report_table = pd.DataFrame()
        # Create pairs table.
        # Used in a precise mode to collect information of insertion coordinates in reference.
        self.pairs_df = self._pairs_table_init()

        # Helper structures:
        # Length of reference contigs
        self.ref_len = dict()
        # Keep track an unique indeces of the reads
        self._index = 0
        # Dictionary of depth attriduted to unclipped reads in given positions
        # unclipped reads do not contain "S" in CIGAR string
        self.unclipped_depth = {}
        # Depth of coverage from clipped reads on the position of junctions.
        # Only clipped reads that DO NOT have junction in the position of interest are count.
        # Required to separate close insertions of IS elements.
        self.cl_read_cov_overlap = {}
        # Boundaries of clipped reads
        self.boundaries = list()
        # Coordinates of provided IS elements. IS name => [chrom, start, stop]
        self.is_coords = dict()
        # List of lengths for matched segments
        # Used to calculate correction coefficients
        self.match_lengths = list()

        # Initialize nested dictionaries as for each dictionary we need contig->position id.
        for contig_i in range(len(self.aln.references)):
            self.ref_len[self.aln.references[contig_i]] = self.aln.lengths[contig_i]
            self.unclipped_depth[self.aln.references[contig_i]] = {}
            self.cl_read_cov_overlap[self.aln.references[contig_i]] = {}

        # Parameters:
        # Show junctions only with this frequency or more
        self.cutoff = 0.005
        # Minimum match in sequences
        self.min_match = 150
        # Average read length
        self.av_read_len = 150
        # Total length of reads
        self.read_lengths = 0
        # Number of analyzed reads
        self.n_reads_analyzed = 0
        # Minimum length of clipped part to use in BLAST
        self.blast_min = 10
        # Maximum expected length of duplication created from the insertion event
        self.max_is_dup_len = 20
        # Minimum identity of clipped read BLAST hit to accept it
        self.blast_min_ident = 75

        # Data folder for circos files
        self.data_folder = './ijump_data/'

    # Initialize a pairs table.
    # Used in a precise mode to collect information of insertion coordinates in reference.
    @staticmethod
    def _pairs_table_init():
        return pd.DataFrame(  # prototype of pairs table
            {'IS_name': ['-'],
             'Position_l': [0],
             'Position_r': [0],
             'Count_mapped_to_IS_l': [0],
             'Count_mapped_to_IS_r': [0],
             'Chrom': ['-']
             }
        )

    # Generate an empty output table
    @staticmethod
    def pairs_table_empty():
        return pd.DataFrame(
            {'Position_l': [],
             'Position_r': [],
             'Count_mapped_to_IS_l': [],
             'Count_mapped_to_IS_r': [],
             'Chrom': [],
             'IS_name': [],
             'Dist': [],
             'N_unclipped_l': [],
             'N_clipped_l': [],
             'N_unclipped_r': [],
             'N_clipped_r': [],
             'N_overlap_l': [],
             'N_overlap_r': [],
             'N_clipped_l_correction': [],
             'N_clipped_r_correction': [],
             'N_overlap_l_correction': [],
             'N_overlap_r_correction': [],
             'N_clipped_l_corrected': [],
             'N_overlap_l_corrected': [],
             'N_clipped_r_corrected': [],
             'N_overlap_r_corrected': [],
             'N_overlap_formula_l': [],
             'N_overlap_formula_r': [],
             'Frequency_l': [],
             'Frequency_r': [],
             'Frequency': [],
             'Depth': []
             }
        )

    # Initialize a report table.
    @staticmethod
    def report_table_init():
        return pd.DataFrame(columns=['IS Name',
                                     'Annotation',
                                     'Chromosome',
                                     'Start',
                                     'Stop',
                                     'Frequency',
                                     'Depth'])

    # Initialize a new copy of clipped reads table.
    @staticmethod
    def _cltbl_init():
        return pd.DataFrame(columns=['ID',
                                     'IS name',
                                     'IS_chrom',
                                     'Read name',
                                     'left pos',  # left position of a clipped segment in a read
                                     'right pos',  # right position of a clipped segment in a read
                                     'clip_position',  # clip position in a reference
                                     'junction_in_read',  # side of a clipped segment connected to a read (l/r)
                                     'reverse',  # is read reverse
                                     'sequence'])  # sequence of a clipped read

    # Initialize a filtered blast output table.
    @staticmethod
    def _blastout_filtered_init():
        return pd.DataFrame(columns=['qseqid',
                                     'sseqid',
                                     'pident',
                                     'length',
                                     'mismatch',
                                     'gapopen',
                                     'qstart',
                                     'qend',
                                     'sstart',
                                     'send',
                                     'evalue',
                                     'bitscore',
                                     'pos_in_ref',
                                     'orientation'])

    # Create summary dataframe.
    def sum_by_reg_tbl_init(self):
        sbrcolumns = ['ann', 'chrom', 'start', 'stop']
        sbrcolumns.extend(list(self.is_coords.keys()))
        return pd.DataFrame(columns=sbrcolumns)

    # Collect information about IS elements.
    def iscollect(self, file):
        logging.info(f'Read file with IS elements: {file}')
        with open(file, 'r') as is_coords_file:
            for coord in is_coords_file.readlines():
                c = coord.split()
                self.is_coords[c[0]] = c[1:]

    # Initialize junction table.
    @staticmethod
    def _jtbl_init(n_rows=0):
        return pd.DataFrame(columns=['ID',
                                     'IS name',
                                     'IS pos',
                                     'IS_chrom',
                                     'Read name',
                                     'Chrom',
                                     'Position',
                                     'Orientation',  # where non-clipped region is relative to position
                                     'Note',
                                     'Locus tag',
                                     'Gene'],
                            index=[i for i in range(n_rows)])

    # For a clipped segment return left, right positions, junction side, coordinate of adjacent non-clipped nucleotide.
    @staticmethod
    def _clboundaries(read):
        positions = read.get_reference_positions(full_length=True)
        clipped_segments = list()
        # Start of clipped segments
        start_clipped_segment = 0
        # Is previous position clipped?
        is_cl_prev = False
        # Is this clipped segment left or right?
        # The segment is "left" if unmapped part of the read is not aligned at current position.
        is_left = None
        # Which part of the read is clipped?
        clipped_part = ''
        # Position of a junction in an aligned part of a read
        junction_pos = 0

        # Collect all clipped segments in a read.
        # Sometimes there are more than one clipped segment (e.g. CIGAR: 30S90M30S).
        for pos_index in range(len(positions)):
            # Check if the position is a start pf a clipped part
            if positions[pos_index] is None and is_cl_prev is False:
                if pos_index == 0:
                    is_left = True
                else:
                    is_left = False
                    clipped_part = 'right'
                    junction_pos = positions[pos_index - 1]
                start_clipped_segment = pos_index + 1
                is_cl_prev = True
            # If current position is (1a) in aligned part of the read or (1b) it is an end of the read
            # and (2) previous position was in unaligned part of the read then (3) assign this position as the end of
            # the clipped part.
            elif (isinstance(positions[pos_index], int) or (pos_index + 1) == len(positions)) and is_cl_prev is True:
                # End of a clipped segment
                end_clipped_segment = pos_index
                if (pos_index + 1) == len(positions):
                    end_clipped_segment = pos_index + 1
                is_cl_prev = False
                if is_left is True:
                    junction_pos = positions[pos_index]
                    clipped_part = 'left'

                clipped_segments.append([start_clipped_segment, end_clipped_segment, clipped_part, junction_pos])
        return clipped_segments

    # return clipped portion of a read
    @staticmethod
    def _clipped_seq(read, left, right):
        return read.query_sequence[(left - 1):right]

    # For each IS element set boundaries where to search clipped reads.
    def set_is_boundaries(self, radius):
        logging.info('Set area near IS elements boundaries to search clipped reads.')
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
                    [self.ref_len[chrom] - radius + start, self.ref_len[chrom], "start", is_name, chrom])
                self.boundaries.append([0, start + radius, "start", is_name, chrom])
            elif start + radius > self.ref_len[chrom]:
                self.boundaries.append([start - radius, self.ref_len[chrom], "start", is_name, chrom])
                self.boundaries.append([0, start + radius - self.ref_len[chrom], "start", is_name, chrom])
            else:
                self.boundaries.append([start - radius, start + radius, "start", is_name, chrom])

            # Check for an end coordinate.
            if stop + radius > self.ref_len[chrom]:
                self.boundaries.append([stop - radius, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append([0, stop + radius - self.ref_len[chrom], "stop", is_name, chrom])
            elif stop - radius < 0:
                self.boundaries.append(
                    [self.ref_len[chrom] - radius + stop, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append([0, stop + radius, "stop", is_name, chrom])
            else:
                self.boundaries.append([stop - radius, stop + radius, "stop", is_name, chrom])

    # Collect clipped reads and check if IS elements are on boundaries of contigs.
    def collect_clipped_reads(self):
        for b in self.boundaries:
            logging.info('Collect clipped reads: ' + ' '.join(str(x) for x in b))

            # Set parameters for clipped reads search
            chrom = b[4]
            start_collection = b[0]
            stop_collection = b[1]
            edge_of_is = b[2]
            name_of_is = b[3]

            # Collect clipped reads
            self._crtable_ungapped(chrom,
                                   start_collection,
                                   stop_collection,
                                   edge_of_is,
                                   name_of_is,
                                   1)

        self.clipped_reads = pd.DataFrame.from_dict(self.clipped_reads_dict, "index")

    # Collect information about coverage that comes from clipped reads outside junction position
    def _cl_read_cov_overlap(self, aln_pairs, chrom):
        if len(aln_pairs) < 3:
            return 0

        read_pos = [a_pair[0] for a_pair in aln_pairs]
        ref_pos = [a_pair[1] for a_pair in aln_pairs]

        for i in read_pos[1:-1]:
            # If nucleotide is not aligned - skip it
            if i is None:
                continue
            else:
                # If the position is junction - skip it
                if ref_pos[i - 1] is None or ref_pos[i + 1] is None:
                    continue
                else:
                    self.cl_read_cov_overlap[chrom][ref_pos[i]] = self.cl_read_cov_overlap[chrom].get(ref_pos[i], 0) + 1

    # Collect clipped reads from the intervals that do not cross boundaries of a contig.
    # The more correct way to collect junction positions would be to find another part of the
    # clipped read in the alignment and take its coordinates. CIGAR strings for both parts could
    # be not mirrored due to some short repeats (1+nt size) near junction positions.
    # direction: 1 => IS->Ref, 0 => Ref->IS (in precise pipeline)
    def _crtable_ungapped(self, chrom, start, stop, edge, is_name, direction):  # generate clipped read table
        # One is added to convert from 0-based to 1-based system
        for read in self.aln.fetch(chrom, start + 1, stop + 1):
            # Add read length to collection of lengths.
            if direction:
                if read.infer_read_length():
                    self.read_lengths += read.infer_read_length()
                    self.n_reads_analyzed += 1

            # Skip unmapped read
            if read.is_unmapped:
                continue

            # Skip if the read is not clipped
            if 'S' not in read.cigarstring:
                continue

            if direction:
                # Collect lengths of read segments that match reference to calculate correction coefficient.
                m_len = [int(x) for x in re.findall(r'(\d+)M', read.cigarstring)]
                # Leave only the longest match from read.
                m_len = max(m_len)
                # Add lengths to collection.
                self.match_lengths.append(m_len)
            else:
                # If it is Ref->IS direction of search:
                # Add coverage from aligned positions of clipped reads that are not junctions.
                self._cl_read_cov_overlap(read.aligned_pairs, read.reference_name)

            # Get clipped segments coordinates from the read
            boundaries = self._clboundaries(read)
            for cl_seg in boundaries:
                # On the IS->Ref search check if read was collected on the correct side of the IS element
                if direction and \
                        not ((cl_seg[2] == "left" and edge == "start") or (cl_seg[2] == "right" and edge == "stop")):
                    continue

                clip_temp = {
                    # Unique read ID
                    'ID': self._index,
                    # IS name
                    'IS name': is_name if direction else '-',
                    # Contig where IS element is located for IS->Ref search and clipped read of Ref->IS
                    'IS_chrom': chrom,
                    'Read name': read.query_name,
                    # Coordinate of clipped segment start
                    'left pos': cl_seg[0],
                    # Coordinate of clipped segment end
                    'right pos': cl_seg[1],
                    # IS clipped segment on "left" from alignment or on "right"
                    'clip_position': cl_seg[2],
                    # Coordinate of junction nucleotide on contig.
                    # At IS side for IS->Ref search and for contig for Ref->IS search
                    'junction_in_read': cl_seg[3],
                    # Is read forward or reverse
                    'reverse': True if read.is_reverse else False,
                    # Sequence of clipped segment
                    'sequence': self._clipped_seq(read, cl_seg[0], cl_seg[1])
                }

                # Add clipped read information to dictionary to build DataFrame.
                # It is faster than append segments-by-segments to the existing DataFrame.
                # As we don't know number of clipped segments we could not create and empty DataFrame of required size.
                if direction:
                    self.clipped_reads_dict[self._index] = clip_temp
                else:
                    self.clipped_reads_bwrd_dict[self._index] = clip_temp

                self._index = self._index + 1

    # Write clipped parts of reads to FASTA file. Use only parts > min_length.
    # direction: 1 => IS->Ref, 0 => Ref->IS
    def _write_cl_fasta(self, cl_fasta_name, min_len, direction):
        if direction:
            cl_table = self.clipped_reads
        else:
            cl_table = self.clipped_reads_bwrd

        with open(cl_fasta_name, 'w') as fasta_file:
            cl_table.index = cl_table['ID']
            for index in cl_table.index:
                if len(cl_table.at[index, 'sequence']) >= min_len:
                    fasta_file.write('>' + str(cl_table.at[index, 'ID']) + '\n')
                    fasta_file.write(str(cl_table.at[index, 'sequence']) + '\n')
                    fasta_file.write('\n')

    # run blast and write output to xml
    # direction: 1 => IS->Ref, 0 => Ref->IS
    def runblast(self, in_file, out_file, direction):
        logging.info('Run BLAST for clipped parts of the reads')
        fasta_file = os.path.join(self.workdir, in_file)
        blast_out_file = os.path.join(self.workdir, out_file)
        self._write_cl_fasta(fasta_file, self.blast_min, direction)

        # Run BLAST
        blastn_cl = NcbiblastnCommandline(query=fasta_file, db=self.ref_name, evalue=0.001, out=blast_out_file,
                                          outfmt=6, word_size=10)
        blastn_cl()

    # Choose left or right coordinate as a clipped junction and orientation relative to junction.
    @staticmethod
    def _choosecoord(qleft, qright, lr):
        qcoord = [qleft, qright]
        qorientation = ['left', 'right']
        coord = int(qcoord[lr == 'left'])
        orientation = qorientation[not (qcoord[1] > qcoord[0]) ^ (lr == 'left')]
        return coord, orientation

    # Parse BALST output.
    # direction: 1=>IS->Ref, 0=>Ref->IS. For Ref->IS we don't need to remove duplicates, we need only one of them
    def parseblast(self, blast_out_file, direction):
        logging.info('Collect information from BLAST')
        blast_out_path = os.path.join(self.workdir, blast_out_file)
        # Check if Blast is not empty
        if os.path.isfile(blast_out_path) and os.path.getsize(blast_out_path) > 0:
            blast_out = pd.read_csv(blast_out_path, sep='\t')
        else:
            logging.info('No BLAST hits were found.')
            raise NoInsertionsFound('No BLAST hits were found.')

        blast_out.columns = ['qseqid',
                             'sseqid',
                             'pident',
                             'length',
                             'mismatch',
                             'gapopen',
                             'qstart',
                             'qend',
                             'sstart',
                             'send',
                             'evalue',
                             'bitscore']

        # Filter only hits with identity [self.blast_min_ident]% or higher. Default: 75%.
        blast_out = blast_out[blast_out['pident'] >= self.blast_min_ident]

        idx_max = blast_out.groupby('qseqid')['bitscore'].transform(max) == blast_out['bitscore']
        # Temporary dataframe for filtering with only best hits by bitscore.
        temp = blast_out[idx_max].copy()

        if direction:
            temp['count'] = temp.groupby('qseqid')['qseqid'].transform('count')
            # Leave only hits with one best hit.
            temp = temp[temp['count'] == 1]
            temp = temp.drop(columns=['count'])
            cl_table = self.clipped_reads
        else:
            temp['rank'] = temp.groupby('qseqid')['bitscore'].rank(method="first", ascending=True)
            # Leave only first best hit.
            temp = temp[temp['rank'] == 1]
            temp = temp.drop(columns=['rank'])
            cl_table = self.clipped_reads_bwrd

        for index in temp.index:
            pos, orient = self._choosecoord(temp.at[index, 'sstart'],
                                            temp.at[index, 'send'],
                                            cl_table.at[temp.at[index, 'qseqid'], 'clip_position'])
            temp.at[index, 'pos_in_ref'] = pos
            temp.at[index, 'orientation'] = orient

        # Check if temp has any entries.
        if temp.size:
            temp['pos_in_ref'] = temp['pos_in_ref'].astype(int)
        else:
            logging.info('No significant BLAST hits.')
            raise NoInsertionsFound('No significant BLAST hits.')

        self.blastout_filtered = temp

    # Check if position close to the IS element
    def _check_is_boundary_proximity(self, chrom, position):
        for b in self.boundaries:
            if b[4] == chrom:
                boundary_width = b[1] - b[0]
                # if b[0] - boundary_width / 2 <= position <= b[1] + boundary_width / 2:  # use doubled boundaries
                if b[0] <= position <= b[1]:
                    return 'IS element', b[3]
        return '-', '-'

    # Cluster positions together to form seed for backwards clipped read search.
    @staticmethod
    def _hclust(X):
        # If only one sample is present – clustering will not work
        if len(X) == 1:
            return [0]
        hcl = AgglomerativeClustering(n_clusters=None, distance_threshold=30, linkage='single'). \
            fit(X.to_numpy().reshape(-1, 1))
        return hcl.labels_

    # Use hierarchical clustering to cluster close positions in the chromosome.
    def make_gene_side_regions(self):
        logging.info('Cluster close positions in the chromosome')

        # Remove positions close to the IS elements boundaries from the analysis
        ref_cl_reads = self.blastout_filtered[['sseqid', 'pos_in_ref']].copy()
        ref_cl_reads = ref_cl_reads.rename(columns={'sseqid': 'Chrom', 'pos_in_ref': 'Position'})
        ref_cl_reads['Note'] = ref_cl_reads.apply(
            lambda x: self._check_is_boundary_proximity(x['Chrom'], x['Position'])[0],
            axis=1
        )
        ref_cl_reads = ref_cl_reads[ref_cl_reads['Note'] == '-']
        ref_cl_reads = ref_cl_reads.drop(columns=['Note'])
        # If no hits point outside IS elements boudaries there is no insertions to find
        if ref_cl_reads.size == 0:
            logging.info("No BLAST hits point oustide IS elements. No significant new insertions could be found.")
            raise NoInsertionsFound(
                "No BLAST hits point oustide IS elements. No significant new insertions could be found."
            )

        ref_cl_reads['Cluster'] = ref_cl_reads. \
            sort_values(by=['Chrom', 'Position']). \
            groupby(['Chrom'])['Position']. \
            transform(self._hclust)

        ref_regions = ref_cl_reads.groupby(['Cluster', 'Chrom']). \
            aggregate(['min', 'max']). \
            reset_index()

        ref_regions.columns = ['Cluster', 'Chrom', 'Position_left', 'Position_right']
        ref_regions = ref_regions.drop(columns=['Cluster'])

        # Extend regions by 5nt if possible
        ref_regions['Position_left'] = ref_regions['Position_left']. \
            apply(lambda x: max(x - 5, 0))
        
        ref_regions['Position_right'] = ref_regions. \
            apply(lambda x: min(x['Position_right'] + 5, self.ref_len[x['Chrom']]), axis=1)

        return ref_regions

    # Collect reads from reference regions (backwards mapping of clipped reads to their IS elements).
    def crtable_bwds(self, ref_regions):
        logging.info('Collect clipped reads from the reference location')
        ref_regions.apply(
            lambda x: self._crtable_ungapped(
                x['Chrom'],
                x['Position_left'],
                x['Position_right'],
                '-',
                '-',
                0
            ),
            axis=1)

        self.clipped_reads_bwrd = pd.DataFrame.from_dict(self.clipped_reads_bwrd_dict, "index")

    # Create table for description of junctions.
    # direction: 1=>IS->Ref, 0=>Ref->IS.
    def call_junctions(self, direction):
        logging.info('Create junction table')
        self.junctions = self._jtbl_init(self.blastout_filtered.shape[0])
        index = 0

        for hit in self.blastout_filtered.itertuples(index=False):
            read_id = hit.qseqid
            if direction:
                pos = hit.pos_in_ref
                chrom = hit.sseqid
                is_name = self.clipped_reads['IS name'][read_id]
                is_chrom = self.clipped_reads['IS_chrom'][read_id]
                is_pos = self.clipped_reads['clip_position'][read_id]
                is_elem_border_mark, _ = self._check_is_boundary_proximity(chrom, pos)
                orientation = hit.orientation
                read_name = self.clipped_reads['Read name'][read_id]
            else:
                pos = self.clipped_reads_bwrd['junction_in_read'][read_id]
                chrom = self.clipped_reads_bwrd['IS_chrom'][read_id]
                is_chrom = hit.sseqid
                is_pos = hit.pos_in_ref
                _, is_name = self._check_is_boundary_proximity(is_chrom, is_pos)
                is_elem_border_mark = '-'
                orientation = self.clipped_reads_bwrd['clip_position'][read_id]
                read_name = self.clipped_reads_bwrd['Read name'][read_id]

            self.junctions.at[index, 'ID'] = read_id
            self.junctions.at[index, 'IS name'] = is_name
            self.junctions.at[index, 'IS_chrom'] = is_chrom
            self.junctions.at[index, 'IS pos'] = is_pos
            self.junctions.at[index, 'Read name'] = read_name
            self.junctions.at[index, 'Chrom'] = chrom
            self.junctions.at[index, 'Position'] = pos
            self.junctions.at[index, 'Orientation'] = orientation
            self.junctions.at[index, 'Locus tag'] = self.gff.gff_pos[chrom][pos][0]
            self.junctions.at[index, 'Gene'] = self.gff.gff_pos[chrom][pos][1]
            self.junctions.at[index, 'Note'] = is_elem_border_mark
            index += 1

        self.junctions = self.junctions.reset_index()

    # Find positions of insertions.
    def search_insert_pos(self):
        logging.info('Serach for junction pairs')
        position_tbl = self.junctions[self.junctions['IS name'] != '-'].copy()
        # It is much better to work with when IS elements collapsed by their copy
        # than to work with each copy separately.
        # Remove copy tags from the IS element names like "_1", "_2".
        position_tbl['IS'] = position_tbl['IS name'].apply(lambda x: re.search(r'(.+)_\d+', x).group(1))
        position_tbl = position_tbl.groupby(['Chrom', 'Position', 'IS', 'Orientation'])['Position']. \
            count(). \
            reset_index(name='Counts')

        # Collect dataframes for pairs of junctions (or orphan junctions) that should mark IS elements insertions.
        is_pairs_collection = []

        # Find pairs
        for chrom in position_tbl['Chrom'].drop_duplicates().tolist():
            # Take IS elements only from the selected chromosome.
            position_tbl_chrom = position_tbl.query('Chrom == @chrom')
            for is_name in position_tbl_chrom['IS'].drop_duplicates().tolist():
                positions_left = position_tbl_chrom.query(
                    'Orientation == "left" & IS == @is_name'
                ).sort_values('Position')
                positions_left_pos = positions_left['Position'].to_numpy()
                positions_left_counts = positions_left['Counts'].to_numpy()
                positions_right = position_tbl_chrom.query(
                    'Orientation == "right" & IS == @is_name'
                ).sort_values('Position')
                positions_right_pos = positions_right['Position'].to_numpy()
                positions_right_counts = positions_right['Counts'].to_numpy()

                logging.info(f'Find pairs for {is_name} and {chrom} contig ')
                # Calculate table of pairs
                pair_tbl_chunk = find_pairs(
                    positions_left_pos,
                    positions_right_pos,
                    positions_left_counts,
                    positions_right_counts,
                    self.ref_len[chrom],
                    self.max_is_dup_len,
                    chrom
                )

                pair_tbl_chunk['IS_name'] = is_name

                is_pairs_collection.append(pair_tbl_chunk)

        # Concatenate all pair tables into one table.
        self.pairs_df = pd.concat(is_pairs_collection, ignore_index=True)

    # Count depth at the position using only unclipped reads.
    # Input is a data frame with columns 'Position' and 'Chrom'.
    def count_depth_unclipped(self, position_tbl):
        logging.info('Count depth attributed to unclipped reads')
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
                elif ('S' not in read.cigarstring) and ('H' not in read.cigarstring):
                    # Test if position is near the end of the read. If it is near skip the read as
                    # it is impossible to distinguish unmapped
                    read_edges = np.array([read.aligned_pairs[0][1], read.aligned_pairs[-1][1]])
                    if np.min(np.abs(pos - read_edges)) > ins_pos_distance:
                        self.unclipped_depth[chrom][pos] = \
                            self.unclipped_depth[chrom].get(pos, 0) + 1

    # Write the full set of output files with headers and zero data rows.
    # A run that finds nothing is a successful run, so it writes the same
    # file set a run that finds something would.
    def _write_empty_outputs(self, mode, outdir):
        self._cltbl_init().to_csv(os.path.join(outdir, "reads.txt"), sep='\t', index=False)
        self._jtbl_init().to_csv(os.path.join(outdir, "ijump_junctions.txt"), sep='\t', index=False)
        self.sum_by_reg_tbl_init().to_csv(os.path.join(outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)
        self.report_table_init().to_csv(os.path.join(outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)

        if mode == EstimationMode.PRECISE:
            self.pairs_table_empty().to_csv(os.path.join(outdir, "ijump_junction_pairs.txt"), sep='\t', index=False)

    # Run the full average/precise pipeline and write every output file it
    # produces along the way. A run that finds nothing to report is not an
    # error: it is signalled through the returned RunResult rather than by
    # letting NoInsertionsFound cross this method's boundary.
    def run(self, mode):
        outdir = os.path.dirname(self.pairs_df_path)

        try:
            # Collect clipped reads information.
            self.collect_clipped_reads()
            # If clipped reads will not be found -> stop the workflow.
            check_data_presence_in_df(self.clipped_reads, 'No clipped reads were found.')
            self.clipped_reads.to_csv(os.path.join(outdir, "reads.txt"), sep='\t', index=False)

            # Run BLAST to search insertion positions in Reference.
            # 1 - search in IS->Reference direction.
            self.runblast('cl.fasta', 'cl_blast.out', 1)
            self.parseblast('cl_blast.out', 1)

            # Read GFF file.
            self.gff.readgff()

            if mode == EstimationMode.AVERAGE:
                # Make a table of observed junction positions
                self.call_junctions(1)

                # Check if any junction is present. If not - stop the workflow.
                check_junctions_presence(self.junctions, outdir, mode)

                # Count reads supporting IS elements insertions for each IS element and each GE
                # Reformat GFF representation
                self.gff.pos_to_ann()
                self.summary_junctions_by_region()
                self.sum_by_region.to_csv(os.path.join(outdir, "ijump_sum_by_reg.txt"), sep='\t', index=False)

                # Make a report table of assessed insertion frequencies in each GE
                self.report_average()
                self.report_table.to_csv(os.path.join(outdir, "ijump_report_by_is_reg.txt"), sep='\t', index=False)
            elif mode == EstimationMode.PRECISE:
                # Make table of regions in the reference genome where extract clipped reads for backwards assignment
                reference_regions = self.make_gene_side_regions()
                reference_regions.to_csv(os.path.join(self.workdir, 'reference_regions.tsv'), sep='\t', index=False)

                # Collect clipped reads at the insertion positions
                # found during forward (IS->Reference) search.
                self.crtable_bwds(reference_regions)

                # Run BLAST to search for positions of found reads.
                check_data_presence_in_df(self.clipped_reads_bwrd, 'No clipped reads were found '
                                                                    'near estimated insertion sites.')
                self.runblast('cl_bwrd.fasta', 'cl_blast_bwds.out', 0)

                # Collect BLAST results.
                self.parseblast('cl_blast_bwds.out', 0)

                # Format results as junction table to attribute reads and their junction positions to the IS elements.
                self.call_junctions(0)

                # Check if any junction is present. If not - stop the workflow.
                check_junctions_presence(self.junctions, outdir, mode)

                # Find pairs of junctions that should indicate insertion positions of both edges of IS element.
                self.search_insert_pos()

                # Filter Junction pairs so at least one of the pair is in the "reference_regions" table.
                self.pairs_df = filter_pairs(self.pairs_df, reference_regions)

                # Check if any pair was produced
                check_data_presence_in_df(self.pairs_df, 'No pairs were found.')

                # Count depth of unclipped reads to have a background depth of coverage
                # Preparation.
                self.pairs_df['Dist'] = self.pairs_df[['Position_l', 'Position_r']]. \
                    apply(
                        lambda clustered_pos: interpos_distance(clustered_pos.Position_l, clustered_pos.Position_r),
                        axis=1
                    )

                positions = pd.concat(
                    [
                        self.pairs_df[['Position_l', 'Chrom', 'Dist']].
                            rename(columns={'Position_l': 'Position'}),
                        self.pairs_df[['Position_r', 'Chrom', 'Dist']].
                            rename(columns={'Position_r': 'Position'})
                    ],
                    axis=0
                ).drop_duplicates()

                # The depth count itself
                self.count_depth_unclipped(positions)
                self.pairs_df.drop(columns='Dist')

                # Make an estimate of insertion frequency
                logging.info('Estimate insertion frequencies.')
                self.pairs_df = frequency_estimation.estimate_frequencies(
                    self.pairs_df, self.clipped_reads_bwrd, self.unclipped_depth,
                    self.cl_read_cov_overlap, self.match_lengths, self.read_lengths, self.n_reads_analyzed,
                )

                # Convert coordinates from 0-base to 1-base
                self.pairs_df['Position_l'] = convert_zero_one_base(self.pairs_df['Position_l'].tolist())
                self.pairs_df['Position_r'] = convert_zero_one_base(self.pairs_df['Position_r'].tolist())
                self.pairs_df.to_csv(self.pairs_df_path, sep='\t', index=False)
        except NoInsertionsFound as e:
            message = str(e)
            logging.info(message)
            self._write_empty_outputs(mode, outdir)
            return RunResult(insertions_found=False, message=message)

        return RunResult(insertions_found=True)

    # make summary table
    def summary_junctions_by_region(self):
        logging.info('Create summary table by region')
        junc_temp = self.junctions.loc[self.junctions['Note'] != 'IS element']
        f_columns = ['ann', 'chrom', 'start', 'stop']
        f_columns.extend(list(self.is_coords.keys()))
        for i in range(len(junc_temp)):
            pos = junc_temp.iloc[i]['Position']
            chrom = junc_temp.iloc[i]['Chrom']
            for ann_id, item in self.gff.ann_pos[chrom].items():  #
                if item[2] <= pos <= item[3]:
                    if ann_id not in self.sum_by_region.index:
                        columns = ['ann_id', 'ann', 'chrom', 'start', 'stop']
                        columns.extend(list(self.is_coords.keys()))
                        temp = self.sum_by_reg_tbl_init()
                        temp.at[0, 'ann_id'] = ann_id
                        temp.at[0, 'ann'] = item[0]
                        temp.at[0, 'chrom'] = item[1]
                        temp.at[0, 'start'] = item[2]
                        temp.at[0, 'stop'] = item[3]
                        temp.at[0, self.is_coords.keys()] = 0
                        temp.at[0, junc_temp.iloc[i]['IS name']] = 1
                        temp = temp.set_index('ann_id')
                        self.sum_by_region = self.sum_by_region.append(temp, sort=True)
                    else:
                        is_name = junc_temp.iloc[i]['IS name']
                        self.sum_by_region.loc[ann_id, is_name] += 1
                    break
        self.sum_by_region = self.sum_by_region[f_columns]

    # Calculate average depth of the region.
    @lru_cache(maxsize=128)
    def _av_depth(self, chrom, start, stop):
        # aln_depth = self.aln.count_coverage(chrom, start, stop)
        # depth = sum(map(sum, aln_depth))
        # return depth / len(aln_depth[0])  # average depth of the region
        c = pysamstats.load_coverage(self.aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000)
        return mean(c.reads_all)

    # Create report by IS and region.
    def report_average(self):
        logging.info("Create report table")
        self.min_match = min(self.match_lengths)  # find minimum match length
        self.av_read_len = self.read_lengths / self.n_reads_analyzed  # find average read length
        self.report_table = pd.melt(
            self.sum_by_region,
            id_vars=('ann', 'chrom', 'start', 'stop'),
            var_name='IS Name',
            value_name='count'
        )

        # Drop zero intervals.
        self.report_table['drop'] = self.report_table.apply(lambda x: 0 if x['stop'] - x['start'] > 0 else 1, axis=1)
        self.report_table = self.report_table[self.report_table['drop'] == 0]
        self.report_table.drop(columns='drop', inplace=True)
        self.report_table.sort_values(by=['start', 'stop'], inplace=True)
        self.report_table = self.report_table[self.report_table['count'] > 0]

        # Add depth.
        self.report_table['Depth'] = self.report_table.apply(
            lambda x: self._av_depth(x['chrom'], x['start'], x['stop']),
            axis=1
        )

        self.report_table['Frequency'] = self.report_table.apply(
            lambda x: round((x['count'] / 2 * (1 + self.blast_min / self.av_read_len)) / (
                    x['Depth'] * (1 - self.min_match / self.av_read_len)), 4),
            axis=1
        )

        self.report_table = self.report_table[['IS Name', 'ann', 'chrom', 'start', 'stop', 'Frequency', 'Depth']]
        self.report_table.columns = ['IS Name', 'Annotation', 'Chromosome', 'Start', 'Stop', 'Frequency', 'Depth']
