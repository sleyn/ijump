# Implement resolve position conflict function when two IS elements are inserted at the same position.

import string
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from Bio.Blast.Applications import NcbiblastnCommandline
import random
import re
import os
import gff
from functools import lru_cache
import pysamstats
from statistics import mean
from sklearn.cluster import AgglomerativeClustering
import logging


# specify class for clipped reads
class ISClipped:
    def __init__(self, aln, ref_name, gff_name, workdir):
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

        # Circos parameters and helper structures:
        # Colors used for IS elements and contig representation
        self._cirocs_colors = ('green',
                               'red',
                               'blue',
                               'purple',
                               'orange',
                               'yellow',
                               'grey')
        # Colors assigned to each chromosome
        self._ref_colours = dict()
        # Colours assigned to each IS element
        self._is_colours = dict()
        # Data folder for circos files
        self.data_folder = './ijump_data/'
        # Id of the session (used in a data folder and config names)
        self.session_id = ''

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
        blast_out = pd.read_csv(os.path.join(self.workdir, blast_out_file), sep='\t')

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

    # Make clusters of left and right insertions junctions from positions.
    # Outputs table of right and left positions of of junctions pairs with counts of clipped
    # reads accounted to each junction.
    # Similar results could be achieved by KNN search, but this algorithm shows slightly better performance on
    # tests.
    @staticmethod
    def _find_pair(pos_l, pos_r, pos_l_count, pos_r_count, chrom_len, max_is_dup_len, chrom):
        # Check if both left and right junctions present. If not - process just present part of junctions.
        if pos_l.size == 0 or pos_r.size == 0:
            n_pairs = pos_l.size + pos_r.size
            pairs_df = pd.DataFrame(
                {'Position_l': [0] * n_pairs,
                 'Position_r': [0] * n_pairs,
                 'Count_mapped_to_IS_l': [0] * n_pairs,
                 'Count_mapped_to_IS_r': [0] * n_pairs,
                 'Chrom': chrom}
            )

            if pos_r.size == 0:
                for pos_l_index, pos in enumerate(pos_l):
                    pairs_df.iloc[pos_l_index, :] = [
                        pos,
                        0,
                        pos_l_count[pos_l_index],
                        0,
                        chrom
                    ]
            else:
                for pos_r_index, pos in enumerate(pos_r):
                    pairs_df.iloc[pos_r_index, :] = [
                        pos,
                        0,
                        pos_r_count[pos_r_index],
                        0,
                        chrom
                    ]

            return pairs_df

        # Store close positions in the matrix where rows are left positions and columns are right positions.
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
            closeness_matrix[pos_index] = (np.ones_like(pos_r) * (np.abs(pos_r - pos) < max_is_dup_len)).astype(np.int0)

        # Assign clusters and sort in each cluster by junction representation in descending order.

        # Build dataframe to populate pairs.
        # We will use maximum number of rows (if all positions do not have pairs).
        n_pairs = np.sum(closeness_matrix.shape)

        pairs_df = pd.DataFrame(
            {'Position_l': [0] * n_pairs,
             'Position_r': [0] * n_pairs,
             'Count_mapped_to_IS_l': [0] * n_pairs,
             'Count_mapped_to_IS_r': [0] * n_pairs,
             'Chrom': chrom}
        )

        # Build clusters of close positions.
        # Clusters are attributed to the left joints.
        cluster_ids = np.zeros(len(closeness_matrix))
        cluster_cur_id = 0
        column_index = 0
        # Itrerate through all closeness_matrix columns or until all left joints will be assigned to clusters.
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
            cluster_pos_l = \
                cluster_pos_l[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            pos_l[np.where(cluster_ids == cluster_id)] = cluster_pos_l

            # Sort closeness sub-matrix.
            cluster_closeness_matrix = closeness_matrix[np.where(cluster_ids == cluster_id)]
            cluster_closeness_matrix = \
                cluster_closeness_matrix[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            closeness_matrix[np.where(cluster_ids == cluster_id), :] = cluster_closeness_matrix

            # Sort counsts sub-list.
            cluster_pos_l_count = pos_l_count[np.where(cluster_ids == cluster_id)]
            cluster_pos_l_count = \
                cluster_pos_l_count[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            pos_l_count[np.where(cluster_ids == cluster_id)] = cluster_pos_l_count

        # Collect right indexes.
        pos_r_orphan = np.arange(pos_r.size)

        # Populate pairs table.
        for pos_l_index, pos_l_cur in enumerate(pos_l):
            if np.sum(closeness_matrix[pos_l_index, :]):
                # Find right index that has minimum difference in counts.
                # Penalize non-cluster items difference by 10000.
                pos_r_index = np.argmin(
                    np.abs(pos_r_count - pos_l_count[pos_l_index]) + ~(closeness_matrix[pos_l_index, :] == 1) * 10000
                )
                pairs_df.iloc[pos_l_index, :] = [
                    pos_l_cur,
                    pos_r[pos_r_index],
                    pos_l_count[pos_l_index],
                    pos_r_count[pos_r_index],
                    chrom
                ]

                closeness_matrix[:, pos_r_index] = 0

                pos_r_orphan[pos_r_index] = -1

            # Write orhphan peaks.
            else:
                pairs_df.iloc[pos_l_index, :] = [
                    pos_l_cur,
                    0,
                    pos_l_count[pos_l_index],
                    0,
                    chrom
                ]

        df_offset = len(pos_l)

        # Add right orphan peaks.
        for shift, pos_r_index_orphan in enumerate(pos_r_orphan[pos_r_orphan != -1]):
            pairs_df.iloc[df_offset + shift, :] = [
                0,
                pos_r[pos_r_index_orphan],
                0,
                pos_r_count[pos_r_index_orphan],
                chrom
            ]

        # Remove empty rows.
        pairs_df = pairs_df.query('Position_l > 0 or Position_r > 0')

        return pairs_df

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

        #
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
                pair_tbl_chunk = self._find_pair(
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

    # Calculate frequency of IS insertion based on frequencies of boundaries junctions.
    @staticmethod
    def _calc_freq_precise(freq_l, freq_r):
        if freq_l == 0:
            return freq_r
        elif freq_r == 0:
            return freq_l
        else:
            return (freq_l + freq_r) / 2

    # Make read count matrices.
    @staticmethod
    def _read_count_mtx(pairs_df, orientation):
        if orientation == 'left':
            pos_df = pairs_df.query('Position_l > 0').copy()
            pos_df = pos_df.rename(columns={'Position_l': 'Position', 'Count_mapped_to_IS_l': 'Count'})
        elif orientation == 'right':
            pos_df = pairs_df.query('Position_r > 0').copy()
            pos_df = pos_df.rename(columns={'Position_r': 'Position', 'Count_mapped_to_IS_r': 'Count'})
        else:
            logging.error('Error: the parameter should be "left" or "right"')
            exit(1)

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
    @staticmethod
    def _resore_orig_counts(counts_mtx, original_rc_counts, pos_dict):
        counts_mtx /= counts_mtx.sum(1).reshape(-1, 1)

        for chr in pos_dict.keys():
            for junct_pos in pos_dict[chr].keys():
                counts_mtx[pos_dict[chr][junct_pos]] *= original_rc_counts[chr].get(junct_pos, 0)

        return counts_mtx

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

    def _add_total_depth(self, depth_l, depth_r):
        if depth_l == 0:
            return depth_r
        elif depth_r == 0:
            return depth_l
        else:
            return (depth_l + depth_r) / 2

    # Add test of excessive reads count that were not mapped to IS elements.
    # The expected difference between N_clipped_l + N_clipped_r (clipped reads collected at junction)
    # and Count_mapped_to_IS_l + Count_mapped_to_IS_r (clipped reads collected at junction that were able
    # to be mapped to IS elements) should be drawn from the uniform distribution and equal
    # P(X <= min_match) = (N_clipped_l + N_clipped_r) * min_match/av_read_len
    # We can compare observed and expected differences with Fisher exact test.
    def fisher_test_clr_number(self, observation):
        sum_clr_mapped_to_is = observation.Count_mapped_to_IS_l + observation.Count_mapped_to_IS_r
        sum_clr_count_at_jnc = observation.N_clipped_l + observation.N_clipped_r

        # We do not expect number of reads mapped to the IS element exceed
        # total number of reads found at the region.
        if sum_clr_count_at_jnc - sum_clr_mapped_to_is < 0:
            return 0

        contingency_table = np.array(
            [
                [sum_clr_mapped_to_is, sum_clr_count_at_jnc * (1 - self.min_match / self.av_read_len)],
                [sum_clr_count_at_jnc - sum_clr_mapped_to_is, sum_clr_count_at_jnc * self.min_match / self.av_read_len]
            ]
        )
        contingency_table.round(0).astype(int)
        _, pvalue = fisher_exact(contingency_table)
        return pvalue

    # Assess frequency of the insertion in population.
    def assess_isel_freq(self):
        # Setup calculation of number of reads supporting each position count
        # for each IS element.
        # 1: Count reads for right and left positions that came directly from positions.
        # Caveat - they do not have information about corresponding IS elements.
        logging.info('Estimate insertion frequencies.')
        original_rc_l_df = self.clipped_reads_bwrd. \
            query('clip_position == "left"'). \
            groupby(['junction_in_read', 'IS_chrom'], as_index=False)['Read name']. \
            count(). \
            rename(columns={'Read name': 'Count', 'IS_chrom': 'Chrom'})

        original_rc_l = {}
        for chrom in self.ref_len.keys():
            original_rc_l[chrom] = {}

        for rc_row_l in original_rc_l_df.itertuples():
            original_rc_l[rc_row_l.Chrom][rc_row_l.junction_in_read] = rc_row_l.Count

        original_rc_r_df = self.clipped_reads_bwrd. \
            query('clip_position == "right"'). \
            groupby(['junction_in_read', 'IS_chrom'], as_index=False)['Read name']. \
            count(). \
            rename(columns={'Read name': 'Count', 'IS_chrom': 'Chrom'})

        original_rc_r = {}
        for chrom in self.ref_len.keys():
            original_rc_r[chrom] = {}

        for rc_row_r in original_rc_r_df.itertuples():
            original_rc_r[rc_row_r.Chrom][rc_row_r.junction_in_read] = rc_row_r.Count

        # 2: Make matrix for left positions.
        pos_l, is_names_l, counts_l = self._read_count_mtx(self.pairs_df, 'left')

        # 3: Make matrix for right positions.
        pos_r, is_names_r, counts_r = self._read_count_mtx(self.pairs_df, 'right')

        # Calculate proportions of reads for each IS for each conflicting position
        # and split reads supporting position from the clipped_reads_bwrd table.
        # 1: left matrix
        counts_l = self._resore_orig_counts(counts_l, original_rc_l, pos_l)

        # 2: right matrix
        counts_r = self._resore_orig_counts(counts_r, original_rc_r, pos_r)

        # Collect numbers of reads at positions for left junctions.
        self.pairs_df['N_unclipped_l'] = self.pairs_df.apply(
            lambda pos: self.unclipped_depth[pos.Chrom].get(pos.Position_l, 0),
            axis=1
        )
        self.pairs_df['N_clipped_l'] = self.pairs_df.apply(
            lambda pair: counts_l[
                pos_l[pair.Chrom][pair.Position_l], is_names_l[pair.IS_name]] if pair.Position_l > 0 else 0,
            axis=1
        )

        # Collect numbers of reads at positions for right junctions.
        self.pairs_df['N_unclipped_r'] = self.pairs_df.apply(
            lambda pos: self.unclipped_depth[pos.Chrom].get(pos.Position_r, 0),
            axis=1
        )
        self.pairs_df['N_clipped_r'] = self.pairs_df.apply(
            lambda pair: counts_r[
                pos_r[pair.Chrom][pair.Position_r], is_names_r[pair.IS_name]
            ] if pair.Position_r > 0 else 0,
            axis=1
        )

        # Add coverage from clipped reads that overlap with junction.
        self.pairs_df['N_overlap_l'] = self.pairs_df[['Position_l', 'Chrom']]. \
            apply(lambda x: self.cl_read_cov_overlap[x.Chrom].get(x.Position_l, 0), axis=1)
        self.pairs_df['N_overlap_r'] = self.pairs_df[['Position_r', 'Chrom']]. \
            apply(lambda x: self.cl_read_cov_overlap[x.Chrom].get(x.Position_r, 0), axis=1)

        # Metrics for corrections and tests
        self.min_match = min(self.match_lengths)
        self.av_read_len = self.read_lengths / self.n_reads_analyzed

        # Add corrections for clipped reads.
        self.pairs_df['N_clipped_l_correction'] = self.pairs_df['N_clipped_l'] / \
                                                  (1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_clipped_l']

        self.pairs_df['N_clipped_r_correction'] = self.pairs_df['N_clipped_r'] / \
                                                  (1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_clipped_r']

        self.pairs_df['N_overlap_l_correction'] = self.pairs_df['N_overlap_l'] / \
                                                  (1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_overlap_l']

        self.pairs_df['N_overlap_r_correction'] = self.pairs_df['N_overlap_r'] / \
                                                  (1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_overlap_r']

        self.pairs_df['N_clipped_l_corrected'] = self.pairs_df['N_clipped_l'] + self.pairs_df['N_clipped_l_correction']
        self.pairs_df['N_overlap_l_corrected'] = self.pairs_df['N_overlap_l'] + self.pairs_df['N_overlap_l_correction']
        self.pairs_df['N_clipped_r_corrected'] = self.pairs_df['N_clipped_r'] + self.pairs_df['N_clipped_r_correction']
        self.pairs_df['N_overlap_r_corrected'] = self.pairs_df['N_overlap_r'] + self.pairs_df['N_overlap_r_correction']

        self.pairs_df['N_overlap_formula_l'] = self.pairs_df[
            ['N_overlap_l_corrected', 'N_clipped_r_corrected', 'Position_l', 'Position_r']
        ]. \
            apply(lambda x: x.N_overlap_l_corrected - x.N_clipped_r_corrected
                            if x.Position_r > x.Position_l > 0
                            else x.N_overlap_l_corrected,
                  axis=1)

        self.pairs_df['N_overlap_formula_r'] = self.pairs_df[
            ['N_overlap_r_corrected', 'N_clipped_l_corrected', 'Position_l', 'Position_r']
        ]. \
            apply(lambda x: x.N_overlap_r_corrected - x.N_clipped_l_corrected
                            if x.Position_r > x.Position_l > 0
                            else x.N_overlap_r_corrected,
                  axis=1)

        # Calculate frequency as average between left and right boundaries if present.
        # If not - just by one boundary.
        # 0.1 pseudocount keeps from div/0 error.
        self.pairs_df['Frequency_l'] = self.pairs_df['N_clipped_l_corrected'] / \
                                       (self.pairs_df['N_unclipped_l'] +
                                        self.pairs_df['N_overlap_formula_l'] +
                                        self.pairs_df['N_clipped_l_corrected'] +
                                        0.1
                                        )
        self.pairs_df['Frequency_r'] = self.pairs_df['N_clipped_r_corrected'] / \
                                       (self.pairs_df['N_unclipped_r'] +
                                        self.pairs_df['N_overlap_formula_r'] +
                                        self.pairs_df['N_clipped_r_corrected'] +
                                        0.1
                                        )

        self.pairs_df['Frequency'] = self.pairs_df[['Frequency_l', 'Frequency_r']]. \
            apply(lambda x: self._calc_freq_precise(x[0], x[1]), axis=1)

        # Add total depth column.
        self.pairs_df['Depth'] = self.pairs_df. \
            apply(
            lambda event:
            self._add_total_depth(
                event.N_unclipped_l + event.N_overlap_formula_l + event.N_clipped_l_corrected,
                event.N_unclipped_r + event.N_overlap_formula_r + event.N_clipped_r_corrected
            ),
            axis=1
        )

    # Generate random string.
    # Was used for Circos.
    @staticmethod
    def _rand_str(n, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(n))

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

    # Create Circos files.
    def create_circos_files(self):
        logging.info('Create CIRCOS files')
        while not os.path.exists(self.data_folder):
            os.makedirs(self.data_folder)

        # Karyotype file
        with open(os.path.join(self.data_folder, 'karyotype.txt'), 'w') as karyotype:
            col_ind = 0
            for contig in self.ref_len.keys():
                karyotype.write('chr - ' + contig + ' ' + contig + ' 0 ' + str(self.ref_len[contig]) + ' ' +
                                self._cirocs_colors[col_ind % len(self._cirocs_colors)] + '\n')
                self._ref_colours[contig] = self._cirocs_colors[col_ind % len(self._cirocs_colors)]
                col_ind += 1

        # Text file
        with open(os.path.join(self.data_folder, 'text.txt'), 'w') as text:
            col_ind = 0
            for is_name in self.is_coords.keys():
                text.write(self.is_coords[is_name][0] + ' ' + self.is_coords[is_name][1] + ' ' +
                           self.is_coords[is_name][1] + ' ' + is_name +
                           ' color=vvd' + self._cirocs_colors[col_ind % len(self._cirocs_colors)] + '\n')
                self._is_colours[is_name] = self._cirocs_colors[col_ind % len(self._cirocs_colors)]
                col_ind += 1

            # List to remove duplicates
            text_regions = list()

            # Add regions information.
            for i in range(len(self.report_table)):
                # Draw only lines with cutoff more than specified.
                if self.report_table.iloc[i]['Frequency'] >= self.cutoff:
                    chrom = self.report_table.iloc[i]['Chromosome']
                    pos = self.report_table.iloc[i]['Start']
                    ann = self.report_table.iloc[i]['Annotation']
                    for a in ann.split('<>'):
                        if a in text_regions:
                            continue
                        text_regions.append(a)
                        text.write(chrom + ' ' + str(pos) + ' ' + str(pos) + ' ' + a + '\n')

        # links
        with open(self.data_folder + 'links.txt', 'w') as links:
            for i in range(len(self.report_table)):
                # Draw only lines with cutoff more than specified
                if self.report_table.iloc[i]['Frequency'] >= self.cutoff:
                    is_name = self.report_table.iloc[i]['IS Name']
                    is_chrom, is_start, is_stop = self.is_coords[is_name]
                    j_chrom = self.report_table.iloc[i]['Chromosome']
                    j_pos = str(self.report_table.iloc[i]['Start'])
                    colour = 'l' + self._is_colours[is_name]
                    links.write(is_chrom + ' ' + is_start + ' ' + is_stop + ' ' + j_chrom + ' ' + j_pos + ' ' +
                                j_pos + ' color=' + colour + '\n')

        # Histogram file
        with open(self.data_folder + 'histogram.txt', 'w') as histogram:
            for i in range(len(self.sum_by_region)):
                # Calculate average depth of the region.
                depth = self._av_depth(self.sum_by_region.iloc[i]['chrom'],
                                       self.sum_by_region.iloc[i]['start'],
                                       self.sum_by_region.iloc[i]['stop'], )

                # Recalculate junction counts to depth.
                h_columns = ['chrom', 'start', 'stop']
                h_columns_is = [x for x in self.is_coords.keys()]

                for h in h_columns_is:
                    if depth > 0:
                        if self.sum_by_region.iloc[i][h] / depth / 2 >= self.cutoff:
                            histogram.write(' '.join(self.sum_by_region.iloc[i][h_columns].apply(str)) + ' ' +
                                            ','.join(self.sum_by_region.iloc[i][h_columns_is].apply(
                                                lambda x: round(((x / depth / 2) * 100), 2)).apply(str)) + '\n')
                            break

        # Depth histogram
        with open(self.data_folder + 'depth.txt', 'w') as depth_hist:
            for contig in self.gff.ann_pos:
                for ann_id, ann in self.gff.ann_pos[contig].items():
                    if ann[3] - ann[2] <= 0:
                        continue
                    depth = self._av_depth(ann[1], ann[2], ann[3])
                    depth_hist.write(' '.join([str(x) for x in ann[1:]]) + ' ' + str(depth) + '\n')

        # Write config.
        config_name = os.path.join(self.data_folder, 'circos.conf')
        with open(config_name, 'w') as config:
            script_folder = os.path.dirname(os.path.realpath(__file__))
            logging.info(script_folder)
            conf_template = open(script_folder + '/circos.conf', 'r')
            conf = conf_template.read()
            conf = re.sub('karyotype = XXX', 'karyotype = ' + self.data_folder + 'karyotype.txt', conf)
            conf = re.sub('XXX		#text', self.data_folder + 'text.txt', conf)
            conf = re.sub('XXX		#links', self.data_folder + 'links.txt', conf)
            conf = re.sub('XXX		#histogram', self.data_folder + 'histogram.txt', conf)
            conf = re.sub('XXX		#depth', self.data_folder + 'depth.txt', conf)

            # Make fill_color string for a histogram.
            hist_colors = ''
            for is_name in self.is_coords.keys():
                hist_colors += self._is_colours[is_name] + ', '
            hist_colors = hist_colors[:-2]

            conf = re.sub('XXX		#stacked_colors', hist_colors, conf)

            config.write(conf)
