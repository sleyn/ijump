# Implement resolve position conflict function when two IS elements are inserted at the same position.

import string
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
import random
import re
import os
import gff
from functools import lru_cache
import pysamstats
from statistics import mean
from sklearn.cluster import AgglomerativeClustering


# specify class for clipped reads
class isclipped:
    def __init__(self, aln, ref_name, gff_name):
        self.aln = aln  # pysam file
        self.clipped_reads = self._cltbl_init()
        self.clipped_reads_bwrd = self._cltbl_init()  # clipped reads for backward run for point pipeline
        self._index = 0  # an unique index of the read
        self.ref_name = ref_name  # file with the reference genome in FASTA format
        self.gff_name = gff_name  # file with GFF annotations
        self.gff = gff.gff(self.gff_name)  # read GFF file
        self.blastout_filtered = self._blastout_filtered_init()  # filered blast output
        self.junctions = self._jtbl_init()  # junction table
        self.ref_len = dict()  # length of reference contigs
        self.unclipped_depth = {}  # dictionary of depth attriduted to unclipped reads
        # (do not contain "S" in CIGAR string)
        for contig_i in range(len(self.aln.references)):
            self.ref_len[self.aln.references[contig_i]] = self.aln.lengths[contig_i]  # Populate lengths of contigs
            self.unclipped_depth[self.aln.references[contig_i]] = {}  # Make nested dictionaries for depth counts

        self.boundaries = list()  # boundaries of clipped reads
        self._cirocs_colors = ('green',
                               'red',
                               'blue',
                               'purple',
                               'orange',
                               'yellow',
                               'grey')

        self._ref_colours = dict()  # colors assigned to each chromosome
        self._is_colours = dict()  # colours assigned to each IS element
        self.is_coords = dict()  # IS name => chrom, start, stop
        self.data_folder = './ijump_data/'  # data folder for circos files
        self.session_id = ''  # id of the session (used in a data folder and config names)
        self.sum_by_region = self._sum_by_reg_tbl()
        self.report_table = self._report_table()
        self.cutoff = 0.005  # show junctions only with this frequency or more
        self.min_match = 150  # minimum match in sequences
        self.av_read_len = 150  # average read length
        self.read_lengths = 0  # total length of reads
        self.n_reads_analyzed = 0  # number of analyzed reads
        self.match_lengths = list()  # list of lengths for matched segments
        self.blast_min = 10  # minimum length of clipped part to use in BLAST
        self.outdir = '.'  # output directory
        self.max_is_dup_len = 20  # maximum expected length of duplication created from the insertion event
        self.pairs_df = pd.DataFrame(  # prototype of pairs table
            {'IS_name': ['-'],
             'Position_l': [0],
             'Position_r': [0],
             'Count_l': [0],
             'Count_r': [0],
             'Chrom': ['-']
             }
        )

    # create report table
    @staticmethod
    def _report_table():
        return pd.DataFrame(columns=['IS Name',
                                     'Annotation',
                                     'Chromosome',
                                     'Start',
                                     'Stop',
                                     'Frequency',
                                     'Depth'])

    # initialize a new copy of clipped reads table
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

    # filered blast output
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

    # collect information about IS elements
    def iscollect(self, file):
        print('Read file ' + file)
        is_coords_file = open(file, 'r')
        for coord in is_coords_file.readlines():
            c = coord.split()
            self.is_coords[c[0]] = c[1:]

        is_coords_file.close()

    # initialize junction table
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

    # for a clipped segment return left, right positions, junction side, coordinate of adjacent non-clipped nucleotide
    @staticmethod
    def _clboundaries(read):
        positions = read.get_reference_positions(full_length=True)
        clipped_segments = list()
        left = 0
        right = 0
        is_cl_prev = False  # is previous position clipped
        is_left = None  # is this clipped segment left or right
        lstr = ''  # string (left/right)
        junction_pos = 0  # position of a junction in an aligned part of a read
        for pos in range(len(positions)):  # collect all clipped segments
            if positions[pos] is None and is_cl_prev is False:
                if pos == 0:
                    is_left = True
                else:
                    is_left = False
                    lstr = 'right'
                    junction_pos = positions[pos - 1]
                left = pos + 1
                is_cl_prev = True
            elif (isinstance(positions[pos], int) or (pos + 1) == len(positions)) and is_cl_prev is True:
                right = pos
                if (pos + 1) == len(positions):
                    right = pos + 1
                is_cl_prev = False
                if is_left is True:
                    junction_pos = positions[pos]
                    lstr = 'left'

                clipped_segments.append([left, right, lstr, junction_pos])
        return clipped_segments

    # create summary dataframe
    def _sum_by_reg_tbl(self):
        sbrcolumns = ['ann', 'chrom', 'start', 'stop']
        sbrcolumns.extend(list(self.is_coords.keys()))
        return pd.DataFrame(columns=sbrcolumns)

    # return clipped portion of a read
    @staticmethod
    def _clipped_seq(read, left, right):
        return read.query_sequence[(left - 1):right]

    # collect clipped reads and check if IS elements are on the boundaries of the contigs
    def crtable(self, radius):
        print('Create table of clipped reads')
        for is_name in self.is_coords.keys():
            chrom = self.is_coords[is_name][0]
            start = int(self.is_coords[is_name][1])
            stop = int(self.is_coords[is_name][2])
            if start - radius < 0:
                self.boundaries.append(
                    [self.ref_len[chrom] - radius + start, self.ref_len[chrom], "start", is_name, chrom])
                self.boundaries.append([0, start + radius, "start", is_name, chrom])
            elif start + radius > self.ref_len[chrom]:
                self.boundaries.append([start - radius, self.ref_len[chrom], "start", is_name, chrom])
                self.boundaries.append([0, start + radius - self.ref_len[chrom], "start", is_name, chrom])
            else:
                self.boundaries.append([start - radius, start + radius, "start", is_name, chrom])

            if stop + radius > self.ref_len[chrom]:
                self.boundaries.append([stop - radius, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append([0, stop + radius - self.ref_len[chrom], "stop", is_name, chrom])
            elif stop - radius < 0:
                self.boundaries.append(
                    [self.ref_len[chrom] - radius + stop, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append([0, stop + radius, "stop", is_name, chrom])
            else:
                self.boundaries.append([stop - radius, stop + radius, "stop", is_name, chrom])

        for b in self.boundaries:
            print(' '.join(str(x) for x in b))
            chrom = b[4]
            start_collection = b[0]
            stop_collection = b[1]
            edge_of_is = b[2]
            name_of_is = b[3]
            self._crtable_ungapped(chrom,
                                   start_collection,
                                   stop_collection,
                                   edge_of_is,
                                   name_of_is,
                                   1)

    # collect clipped reads from the intervals that do not cross boundaries of a contig
    # The more correct way to collect junction positions would be to find another part of the
    # clipped read in the alignment and take its coordinates. CIGAR strings for both parts could
    # be not mirrored due to some short repeats (1+nt size) near junction positions.
    # run: 1 => IS->Ref, 0 => Ref->IS (in point pipeline)
    def _crtable_ungapped(self, chrom, start, stop, edge, is_name, run):  # generate clipped read table
        # one is added to convert from 0-based to 1-based system
        for read in self.aln.fetch(chrom, start + 1, stop + 1):
            if run:
                if read.infer_read_length():
                    self.read_lengths += read.infer_read_length()  # add read length to collection of lengths
                    self.n_reads_analyzed += 1

            if read.is_unmapped:
                continue  # skip unmapped read
            elif 'S' not in read.cigarstring:
                continue  # skip not clipped read

            if run:
                m_len = [int(x) for x in re.findall('(\d+)M', read.cigarstring)]  # collect matched lengths
                m_len = max(m_len)  # leave only longest match from read
                self.match_lengths.append(m_len)  # add lengths to collection

            #           for i in range(read.cigarstring.count('S')):
            boundaries = self._clboundaries(read)
            for cl_seg in boundaries:
                # On the IS->Ref search check if read was collected on the correct side of the IS element
                if run and not ((cl_seg[2] == "left" and edge == "start") or (cl_seg[2] == "right" and edge == "stop")):
                    continue

                clip_temp = self._cltbl_init()
                clip_temp.at[0, 'ID'] = self._index

                if run:
                    clip_temp.at[0, 'IS name'] = is_name
                else:
                    clip_temp.at[0, 'IS name'] = '-'

                clip_temp.at[0, 'IS_chrom'] = chrom
                clip_temp.at[0, 'Read name'] = read.query_name
                clip_temp.at[0, 'left pos'] = cl_seg[0]
                clip_temp.at[0, 'right pos'] = cl_seg[1]
                clip_temp.at[0, 'clip_position'] = cl_seg[2]
                clip_temp.at[0, 'junction_in_read'] = cl_seg[3]

                if read.is_reverse:
                    clip_temp.at[0, 'reverse'] = True
                else:
                    clip_temp.at[0, 'reverse'] = False

                clip_temp.at[0, 'sequence'] = self._clipped_seq(read, cl_seg[0], cl_seg[1])

                if run:
                    self.clipped_reads = self.clipped_reads.append(clip_temp)
                else:
                    self.clipped_reads_bwrd = self.clipped_reads_bwrd.append(clip_temp)

                self._index = self._index + 1

    # write clipped reads to fasta file with minimum length
    # run: 1 => IS->Ref, 0 => Ref->IS
    def _write_fasta(self, cl_fasta_name, min_len, run):
        if run:
            cl_table = self.clipped_reads
        else:
            cl_table = self.clipped_reads_bwrd

        fasta_file = open(cl_fasta_name, 'w')
        cl_table.index = cl_table['ID']
        for index in cl_table.index:
            if len(cl_table.at[index, 'sequence']) >= min_len:
                fasta_file.write('>' + str(cl_table.at[index, 'ID']) + '\n')
                fasta_file.write(str(cl_table.at[index, 'sequence']) + '\n')
                fasta_file.write('\n')
        fasta_file.close()

    # run blast and write output to xml
    # run: 1 => IS->Ref, 0 => Ref->IS
    def runblast(self, in_file, out_file, run):
        print('Run BLAST for clipped parts of the reads')
        fasta_file = os.path.join(self.outdir, in_file)
        blast_out_file = os.path.join(self.outdir, out_file)
        self._write_fasta(fasta_file, self.blast_min, run)
        blastn_cl = NcbiblastnCommandline(query=fasta_file, db=self.ref_name, evalue=0.001, out=blast_out_file,
                                          outfmt=6, word_size=10)
        blastn_cl()

    # choose left or right coordinate as a clipped junction and orientation relative to junction
    @staticmethod
    def _choosecoord(qleft, qright, lr):
        qcoord = [qleft, qright]
        qorientation = ['left', 'right']
        coord = int(qcoord[lr == 'left'])
        orientation = qorientation[not (qcoord[1] > qcoord[0]) ^ (lr == 'left')]
        return coord, orientation

    # parse BALST output
    # run: 1=>IS->Ref, 0=>Ref->IS. For Ref->IS we don't need to remove duplicates, we need only one of them
    def parseblast(self, run, blast_out_file):
        print('Collect information from BLAST')
        blast_out = pd.read_csv(os.path.join(self.outdir, blast_out_file), sep='\t')

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

        blast_out = blast_out[blast_out['pident'] >= 75]  # filter only hits with 75% or higher

        idx_max = blast_out.groupby('qseqid')['bitscore'].transform(max) == blast_out['bitscore']
        temp = blast_out[idx_max].copy()  # temporary dataframe for filtering with only best hits by bitscore

        if run:
            temp['count'] = temp.groupby('qseqid')['qseqid'].transform('count')
            temp = temp[temp['count'] == 1]  # leave only hits with one best hit
            temp = temp.drop(columns=['count'])
        else:
            temp['rank'] = temp.groupby('qseqid')['bitscore'].rank(method="first", ascending=True)
            temp = temp[temp['rank'] == 1]  # leave only first best hit
            temp = temp.drop(columns=['rank'])

        if run:
            cl_table = self.clipped_reads
        else:
            cl_table = self.clipped_reads_bwrd

        for index in temp.index:
            pos, orient = self._choosecoord(temp.at[index, 'sstart'],
                                            temp.at[index, 'send'],
                                            cl_table.at[temp.at[index, 'qseqid'], 'clip_position'])
            temp.at[index, 'pos_in_ref'] = pos
            temp.at[index, 'orientation'] = orient

        # if temp has any entries
        if temp.size:
            temp['pos_in_ref'] = temp['pos_in_ref'].astype(int)
        self.blastout_filtered = temp

    # Check if position close to the IS element
    def _check_is_boundary_proximity(self, chrom, position):
        for b in self.boundaries:
            if b[4] == chrom:
                boundary_width = b[1] - b[0]
                if b[0] - boundary_width / 2 <= position <= b[1] + boundary_width / 2:  # use doubled boundaries
                    return 'IS element', b[3]
        return '-', '-'

    # use hierarchical clustering to cluster close positions in the chromosome
    def make_gene_side_regions(self):
        # Remove positions close to the IS elements boundaries from the analysis
        ref_cl_reads = self.blastout_filtered[['sseqid', 'pos_in_ref']]
        ref_cl_reads = ref_cl_reads.rename(columns={'sseqid': 'Chrom', 'pos_in_ref': 'Position'})
        ref_cl_reads['Note'] = ref_cl_reads.apply(
            lambda x: self._check_is_boundary_proximity(x['Chrom'], x['Position'])[0],
            axis=1
        )
        ref_cl_reads = ref_cl_reads[ref_cl_reads['Note'] == '-']
        ref_cl_reads = ref_cl_reads.drop(columns=['Note'])

        # Cluster positions together to form seed for backwards clipped read search
        def hclust(X):
            # If only one sample is present – clustering will not work
            if len(X) == 1:
                return [0]
            hcl = AgglomerativeClustering(n_clusters=None, distance_threshold=30, linkage='single'). \
                fit(X.to_numpy().reshape(-1, 1))
            return hcl.labels_

        ref_cl_reads['Cluster'] = ref_cl_reads.groupby(['Chrom'])['Position']. \
            transform(hclust)

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

    # collect reads from reference regions (backwards mapping of clipped reads to their IS elements)
    def crtable_bwds(self, ref_regions):
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

    # create table for junctions description
    def call_junctions(self, run):
        print('Create junction table')
        self.junctions = self._jtbl_init(self.blastout_filtered.shape[0])
        index = 0

        for hit in self.blastout_filtered.itertuples(index=False):
            read_id = hit.qseqid
            if run:
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
                 'Count_l': [0] * n_pairs,
                 'Count_r': [0] * n_pairs,
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

        # Store close positions in the matrix with rows as left positions and columns are right positions
        # The value is 1 if two positions are closer then max_is_dup_len value
        closeness_matrix = np.zeros((pos_l.size, pos_r.size))

        # Check if any position close to the contig ends
        if pos_r[-1] - pos_l[0] > chrom_len / 2:
            if chrom_len - (pos_r[-1] - pos_l[0]) <= max_is_dup_len:
                closeness_matrix[0, -1] = 1

        if pos_l[-1] - pos_r[0] > chrom_len / 2:
            if chrom_len - (pos_l[-1] - pos_r[0]) <= max_is_dup_len:
                closeness_matrix[-1, 0] = 1

        # Populate closeness matrix
        for pos_index, pos in enumerate(pos_l):
            closeness_matrix[pos_index] = (np.ones_like(pos_r) * (np.abs(pos_r - pos) < max_is_dup_len)).astype(np.int0)

        # Assign clusters and sort in each cluster by junction representation in descending order

        # Build dataframe to populate pairs
        # We will use maximum number of rows (if all positions do not have pairs)
        n_pairs = np.sum([
            closeness_matrix[closeness_matrix.any(1)][:, closeness_matrix.any(0)].shape[0],
            closeness_matrix[closeness_matrix.any(1)][:, closeness_matrix.any(0)].shape[1]
        ])

        pairs_df = pd.DataFrame(
            {'Position_l': [0] * n_pairs,
             'Position_r': [0] * n_pairs,
             'Count_l': [0] * n_pairs,
             'Count_r': [0] * n_pairs,
             'Chrom': chrom}
        )

        # Build clusters of close positions
        # Clusters are attributed to the left joints
        cluster_ids = np.zeros(len(closeness_matrix))
        cluster_cur_id = 0
        column_index = 0
        # Itrerate through all closeness_matrix columns or before all left joints will be assigned to clusters
        while not (column_index >= closeness_matrix.shape[1] or np.all(cluster_ids > 0)):
            # Check if the column not zero (orphan right position)
            if closeness_matrix[:, column_index].any():
                # If any left position has several right positions in proximity
                # unite clusters
                if np.any(closeness_matrix[:, column_index][cluster_ids > 0] == 1):
                    cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id
                else:
                    # If cluster is first or clusters do not overlap add cluster id
                    cluster_cur_id += 1
                    cluster_ids[closeness_matrix[:, column_index] == 1] = cluster_cur_id

            column_index += 1

        # Sort each cluster
        for cluster_id in np.unique(cluster_ids[cluster_ids > 0]):
            # Sort positions sub-list
            cluster_pos_l = pos_l[np.where(cluster_ids == cluster_id)]
            cluster_pos_l = \
                cluster_pos_l[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            pos_l[np.where(cluster_ids == cluster_id)] = cluster_pos_l

            # Sort closeness sub-matrix
            cluster_closeness_matrix = closeness_matrix[np.where(cluster_ids == cluster_id)]
            cluster_closeness_matrix = \
                cluster_closeness_matrix[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            closeness_matrix[np.where(cluster_ids == cluster_id), :] = cluster_closeness_matrix

            # Sort counsts sub-list
            cluster_pos_l_count = pos_l_count[np.where(cluster_ids == cluster_id)]
            cluster_pos_l_count = \
                cluster_pos_l_count[np.argsort(pos_l_count[np.where(cluster_ids == cluster_id)])[::-1]]
            pos_l_count[np.where(cluster_ids == cluster_id)] = cluster_pos_l_count

        # Collect right indexes
        pos_r_orphan = np.arange(pos_r.size)

        # Populate pairs table
        for pos_l_index, pos_l_cur in enumerate(pos_l):
            if np.sum(closeness_matrix[pos_l_index, :]):
                # Find right index that has minimum difference in counts
                # Penalize non-cluster items difference by 10000
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

            # Write orhphan peaks
            else:
                pairs_df.iloc[pos_l_index, :] = [
                    pos_l_cur,
                    0,
                    pos_l_count[pos_l_index],
                    0,
                    chrom
                ]

        df_offset = pos_l_index + 1

        # Add right orphan peaks
        for shift, pos_r_index_orphan in enumerate(pos_r_orphan[pos_r_orphan != -1]):
            pairs_df.iloc[df_offset + shift, :] = [
                0,
                pos_r[pos_r_index_orphan],
                0,
                pos_r_count[pos_r_index_orphan],
                chrom
            ]

        # Remove empty rows
        pairs_df = pairs_df.query('Position_l > 0 or Position_r > 0')

        return pairs_df

    # Find positions of insertions
    def search_insert_pos(self):
        print('Serach for junction pairs')
        position_tbl = self.junctions[self.junctions['IS name'] != '-'].copy()
        # It is much better to work with when IS elements collapsed by their copy then to work with each copy separately
        # remove copy tags from the IS element names like "_1", "_2".
        position_tbl['IS'] = position_tbl['IS name'].apply(lambda x: re.search(r'(.+)_\d+', x).group(1))
        position_tbl = position_tbl.groupby(['Chrom', 'Position', 'IS', 'Orientation'])['Position']. \
            count(). \
            reset_index(name='Counts')

        # collect dataframes for pairs of junctions (or orphan junctions) that should mark IS elements insertions
        is_pairs_collection = []

        for chrom in position_tbl['Chrom'].drop_duplicates().tolist():
            # Take IS elements only from the selected chromsome
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

                print(f'Find pairs for {is_name}')
                # Calculate table
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

        # Concatenate all pair tables into one table
        self.pairs_df = pd.concat(is_pairs_collection, ignore_index=True)

    # Calculate frequency
    @staticmethod
    def _calc_freq_point(freq_l, freq_r):
        if freq_l == 0:
            return freq_r
        elif freq_r == 0:
            return freq_l
        else:
            return (freq_l + freq_r) / 2

    # Make read count matrices
    @staticmethod
    def _read_count_mtx(pairs_df, orientation):
        if orientation == 'left':
            pos_df = pairs_df.query('Position_l > 0').copy()
            pos_df = pos_df.rename(columns={'Position_l': 'Position', 'Count_l': 'Count'})
        elif orientation == 'right':
            pos_df = pairs_df.query('Position_r > 0').copy()
            pos_df = pos_df.rename(columns={'Position_r': 'Position', 'Count_r': 'Count'})
        else:
            print('Error: the parameter should be "left" or "right"')
            exit(1)

        # Dictionary to translate positions (Contig name/Coordinate) to matrix row indeces
        pos = {}
        for chrom in pos_df.Chrom.unique():
            pos[chrom] = {}

        i = 0

        for pos_row in pos_df[['Position', 'Chrom']].drop_duplicates().itertuples():
            pos[pos_row.Chrom][pos_row.Position] = i
            i += 1

        # Dictionary to translate IS names to matrix row indeces
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

    # Translate matrix of read count proporions to the original read counts
    # (calculated from the second pass of clipped reads collection)
    @staticmethod
    def _resore_orig_counts(counts_mtx, original_rc_counts, pos_dict):
        counts_mtx /= counts_mtx.sum(1).reshape(-1, 1)

        for chr in pos_dict.keys():
            for junct_pos in pos_dict[chr].keys():
                counts_mtx[pos_dict[chr][junct_pos]] *= original_rc_counts[chr].get(junct_pos, 0)

        return counts_mtx

    # Count depth at the position using only unclipped reads
    # Input is a data frame with columns 'Position' and 'Chrom'
    def count_depth(self, position_tbl):
        for position in position_tbl.itertuples():
            chrom = position.Chrom
            pos = position.Position
            if pos == 0:
                continue

            for read in self.aln.fetch(chrom, pos, pos + 1):
                if read.is_unmapped:
                    continue  # skip unmapped read
                elif 'S' not in read.cigarstring:
                    self.unclipped_depth[chrom][pos] = \
                        self.unclipped_depth[chrom].get(pos, 0) + 1

    # Assess frequency of the insertion in population
    def assess_isel_freq(self):
        # Setup calculation of number of reads supporting each position count
        # for each IS element
        # 1: Count reads for right and left positions that came directly from positions.
        # Caveat - they do not have information about corresponding IS elements.
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

        # 2: Make matrix for left positions
        pos_l, is_names_l, counts_l = self._read_count_mtx(self.pairs_df, 'left')

        # 3: Make matrix for right positions
        pos_r, is_names_r, counts_r = self._read_count_mtx(self.pairs_df, 'right')

        # Calculate proportions of reads for each IS for each conflicting position
        # and split reads supporting position from the clipped_reads_bwrd table
        # 1: left matrix
        counts_l = self._resore_orig_counts(counts_l, original_rc_l, pos_l)

        # 2: right matrix
        counts_r = self._resore_orig_counts(counts_r, original_rc_r, pos_r)

        # Collect numbers of reads at positions for left junctions
        self.pairs_df['N_unclipped_l'] = self.pairs_df.apply(
            lambda pos: self.unclipped_depth[pos.Chrom].get(pos.Position_l, 0),
            axis=1
        )
        self.pairs_df['N_clipped_l'] = self.pairs_df.apply(
            lambda pair: counts_l[
                pos_l[pair.Chrom][pair.Position_l], is_names_l[pair.IS_name]] if pair.Position_l > 0 else 0,
            axis=1
        )

        # Collect numbers of reads at positions for right junctions
        self.pairs_df['N_unclipped_r'] = self.pairs_df.apply(
            lambda pos: self.unclipped_depth[pos.Chrom].get(pos.Position_r, 0),
            axis=1
        )
        self.pairs_df['N_clipped_r'] = self.pairs_df.apply(
            lambda pair: counts_r[
                pos_r[pair.Chrom][pair.Position_r], is_names_l[pair.IS_name]] if pair.Position_r > 0 else 0,
            axis=1
        )

        # Add corrections for clipped reads
        self.min_match = min(self.match_lengths)
        self.av_read_len = self.read_lengths / self.n_reads_analyzed

        self.pairs_df['N_clipped_l_correction'] = self.pairs_df['N_clipped_l'] * \
                                                  (1 + self.blast_min / self.av_read_len) / (
                                                              1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_clipped_l']
        self.pairs_df['N_clipped_r_correction'] = self.pairs_df['N_clipped_r'] * \
                                                  (1 + self.blast_min / self.av_read_len) / (
                                                              1 - self.min_match / self.av_read_len) - \
                                                  self.pairs_df['N_clipped_r']

        self.pairs_df['N_clipped_l'] = self.pairs_df['N_clipped_l'] + self.pairs_df['N_clipped_l_correction']
        self.pairs_df['N_unclipped_l'] = self.pairs_df['N_unclipped_l']
        self.pairs_df['N_clipped_r'] = self.pairs_df['N_clipped_r'] + self.pairs_df['N_clipped_r_correction']
        self.pairs_df['N_unclipped_r'] = self.pairs_df['N_unclipped_r']

        # Calculate frequency as average between left and right boundaries if present
        # If not - just by one boundary
        # 0.1 pseudocount keeps from div/0 error
        self.pairs_df['Frequency_l'] = self.pairs_df['N_clipped_l'] / \
                                       (self.pairs_df['N_unclipped_l'] + self.pairs_df['N_clipped_l'] + 0.1)
        self.pairs_df['Frequency_r'] = self.pairs_df['N_clipped_r'] / \
                                       (self.pairs_df['N_unclipped_r'] + self.pairs_df['N_clipped_r'] + 0.1)

        self.pairs_df['Frequency'] = self.pairs_df[['Frequency_l', 'Frequency_r']]. \
            apply(lambda x: self._calc_freq_point(x[0], x[1]), axis=1)

    # generate random string
    @staticmethod
    def _rand_str(n, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(n))

    # make summary table
    def summary_junctions_by_region(self):
        print('Create summary table by region')
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
                        temp = self._sum_by_reg_tbl()
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

    # calculate average depth of the region
    @lru_cache(maxsize=128)
    def _av_depth(self, chrom, start, stop):
        # aln_depth = self.aln.count_coverage(chrom, start, stop)
        # depth = sum(map(sum, aln_depth))
        # return depth / len(aln_depth[0])  # average depth of the region
        c = pysamstats.load_coverage(self.aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000)
        return mean(c.reads_all)

    # create report by IS and region
    def report(self):
        print("Create report table")
        self.min_match = min(self.match_lengths)  # find minimum match length
        self.av_read_len = self.read_lengths / self.n_reads_analyzed  # find average read length
        self.report_table = pd.melt(
            self.sum_by_region,
            id_vars=('ann', 'chrom', 'start', 'stop'),
            var_name='IS Name',
            value_name='count'
        )

        # Drop zero intervals
        self.report_table['drop'] = self.report_table.apply(lambda x: 0 if x['stop'] - x['start'] > 0 else 1, axis=1)
        self.report_table = self.report_table[self.report_table['drop'] == 0]
        self.report_table.drop(columns='drop', inplace=True)
        self.report_table.sort_values(by=['start', 'stop'], inplace=True)
        self.report_table = self.report_table[self.report_table['count'] > 0]

        # Add depth
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

    # create Circos files
    def create_circos_files(self):
        print('Create CIRCOS files')
        while not os.path.exists(self.data_folder):
            # self.session_id = self._rand_str(5)
            # self.data_folder = "./data" + self.session_id + "/"
            os.makedirs(self.data_folder)

        # karyotype file
        karyotype = open(self.data_folder + 'karyotype.txt', 'w')
        col_ind = 0
        for contig in self.ref_len.keys():
            karyotype.write('chr - ' + contig + ' ' + contig + ' 0 ' + str(self.ref_len[contig]) + ' ' +
                            self._cirocs_colors[col_ind % len(self._cirocs_colors)] + '\n')
            self._ref_colours[contig] = self._cirocs_colors[col_ind % len(self._cirocs_colors)]
            col_ind += 1
        karyotype.close()

        # text
        text = open(self.data_folder + 'text.txt', 'w')

        col_ind = 0
        for is_name in self.is_coords.keys():
            text.write(self.is_coords[is_name][0] + ' ' + self.is_coords[is_name][1] + ' ' +
                       self.is_coords[is_name][1] + ' ' + is_name +
                       ' color=vvd' + self._cirocs_colors[col_ind % len(self._cirocs_colors)] + '\n')
            self._is_colours[is_name] = self._cirocs_colors[col_ind % len(self._cirocs_colors)]
            col_ind += 1

        text_regions = list()  # to remove duplicates

        # add regions information
        for i in range(len(self.report_table)):
            if self.report_table.iloc[i]['Frequency'] >= self.cutoff:  # draw only lines with cutoff more than specified
                chrom = self.report_table.iloc[i]['Chromosome']
                pos = self.report_table.iloc[i]['Start']
                ann = self.report_table.iloc[i]['Annotation']
                for a in ann.split('<>'):
                    if a in text_regions:
                        continue
                    text_regions.append(a)
                    text.write(chrom + ' ' + str(pos) + ' ' + str(pos) + ' ' + a + '\n')

        text.close()

        # links
        links = open(self.data_folder + 'links.txt', 'w')

        for i in range(len(self.report_table)):
            if self.report_table.iloc[i]['Frequency'] >= self.cutoff:  # draw only lines with cutoff more than specified
                is_name = self.report_table.iloc[i]['IS Name']
                is_chrom, is_start, is_stop = self.is_coords[is_name]
                j_chrom = self.report_table.iloc[i]['Chromosome']
                j_pos = str(self.report_table.iloc[i]['Start'])
                colour = 'l' + self._is_colours[is_name]
                links.write(is_chrom + ' ' + is_start + ' ' + is_stop + ' ' + j_chrom + ' ' + j_pos + ' ' +
                            j_pos + ' color=' + colour + '\n')

        links.close()

        # histogram
        histogram = open(self.data_folder + 'histogram.txt', 'w')

        for i in range(len(self.sum_by_region)):
            depth = self._av_depth(self.sum_by_region.iloc[i]['chrom'],
                                   self.sum_by_region.iloc[i]['start'],
                                   self.sum_by_region.iloc[i]['stop'], )  # average depth of the region

            # recalculate junction counts to depth
            h_columns = ['chrom', 'start', 'stop']
            h_columns_is = [x for x in self.is_coords.keys()]

            for h in h_columns_is:
                if depth > 0:
                    if self.sum_by_region.iloc[i][h] / depth / 2 >= self.cutoff:
                        histogram.write(' '.join(self.sum_by_region.iloc[i][h_columns].apply(str)) + ' ' +
                                        ','.join(self.sum_by_region.iloc[i][h_columns_is].apply(
                                            lambda x: round(((x / depth / 2) * 100), 2)).apply(str)) + '\n')
                        break

        histogram.close()

        # depth histogram
        depth_hist = open(self.data_folder + 'depth.txt', 'w')

        for contig in self.gff.ann_pos:
            for ann_id, ann in self.gff.ann_pos[contig].items():
                if ann[3] - ann[2] <= 0:
                    continue
                depth = self._av_depth(ann[1], ann[2], ann[3])
                depth_hist.write(' '.join([str(x) for x in ann[1:]]) + ' ' + str(depth) + '\n')

        depth_hist.close()

        # write config
        config_name = self.data_folder + 'circos.conf'
        config = open(config_name, 'w')
        script_folder = os.path.dirname(os.path.realpath(__file__))
        print(script_folder)
        conf_template = open(script_folder + '/circos.conf', 'r')
        conf = conf_template.read()
        conf = re.sub('karyotype = XXX', 'karyotype = ' + self.data_folder + 'karyotype.txt', conf)
        conf = re.sub('XXX		#text', self.data_folder + 'text.txt', conf)
        conf = re.sub('XXX		#links', self.data_folder + 'links.txt', conf)
        conf = re.sub('XXX		#histogram', self.data_folder + 'histogram.txt', conf)
        conf = re.sub('XXX		#depth', self.data_folder + 'depth.txt', conf)

        # make fill_color string for a histogram
        hist_colors = ''
        for is_name in self.is_coords.keys():
            hist_colors += self._is_colours[is_name] + ', '
        hist_colors = hist_colors[:-2]

        conf = re.sub('XXX		#stacked_colors', hist_colors, conf)

        config.write(conf)
