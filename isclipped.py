import string
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
import random
import re
import os
import gff
from functools import lru_cache
import pysamstats
from statistics import mean

# specify class for clipped reads

class isclipped:
    def __init__(self, aln, ref_name, gff_name):
        self.aln = aln  # pysam file
        self.clipped_reads = self._cltbl_init()
        self._index = 0  # an unique index of the read
        self.ref_name = ref_name  # file with the reference genome in FASTA format
        self.gff_name = gff_name  # file with GFF annotations
        self.gff = gff.gff(self.gff_name)       # read GFF file
        self.blastout_filtered = pd.DataFrame(columns=['qseqid',  # filered blast output
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
                                                       'pos in ref',
                                                       'orientation'])
        self.junctions = self._jtbl_init()
        self.ref_len = dict()
        for contig_i in range(len(self.aln.references)):
            self.ref_len[self.aln.references[contig_i]] = self.aln.lengths[contig_i]

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
        self.cutoff = 0.005     # show junctions only with this frequency or more
        self.min_match = 150     # minimum match in sequences
        self.av_read_len = 150  # average read length
        self.read_lengths = 0 # total length of reads
        self.n_reads_analyzed = 0 # number of analyzed reads
        self.match_lengths = list() # list of lengths for matched segments
        self.blast_min = 10 # minimum length of clipped part to use in BLAST
        self.outdir = '.' # output directory

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

    # collect information about IS elements
    def iscollect(self, file):
        print('Read file ' + file)
        is_coords_file = open(file, 'r')
        for coord in is_coords_file.readlines():
            c = coord.split()
            self.is_coords[c[0]] = c[1:]

        is_coords_file.close()

    # initialize a new copy of clipped reads table
    @staticmethod
    def _cltbl_init():
        return pd.DataFrame(columns=['ID',
                                     'IS name',
                                     'IS chrom',
                                     'Read name',
                                     'left pos',  # left position of a clipped segment in a read
                                     'right pos',  # right position of a clipped segment in a read
                                     'clip position',  # clip position in a reference
                                     'junction in read',  # side of a clipped segment connected to a read (l/r)
                                     'reverse',  # is read reverse
                                     'sequence'])  # sequence of a clipped read

    # initialize junction table
    @staticmethod
    def _jtbl_init():
        return pd.DataFrame(columns=['ID',
                                     'IS name',
                                     'IS pos',
                                     'IS chrom',
                                     'Read name',
                                     'Chrom',
                                     'Position',
                                     'Orientation',  # where non-clipped region is relative to position
                                     'Note',
                                     'Locus tag',
                                     'Gene'])

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
        for pos in range(len(positions)):               # collect all clipped segments
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
                self.boundaries.append([self.ref_len[chrom] - radius + start, self.ref_len[chrom], "start", is_name, chrom])
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
                self.boundaries.append([self.ref_len[chrom] - radius + stop, self.ref_len[chrom], "stop", is_name, chrom])
                self.boundaries.append([0, stop + radius, "stop", is_name, chrom])
            else:
                self.boundaries.append([stop - radius, stop + radius, "stop", is_name, chrom])

        for b in self.boundaries:
            print(' '.join(str(x) for x in b))
            self._crtable_ungapped(b[4], b[0], b[1], b[2], b[3])

    # collect clipped reads from the intervals that do not cross boundaries of a contig
    def _crtable_ungapped(self, chrom, start, stop, edge, is_name):  # generate clipped read table
        for read in self.aln.fetch(chrom, start + 1,
                                   stop + 1):  # one is added to convert from 0-based to 1-based system
            if read.infer_read_length():
                self.read_lengths += read.infer_read_length()      # add read length to collection of lengths
                self.n_reads_analyzed += 1

            if read.is_unmapped:
                continue            # skip unmapped read
            elif 'S' not in read.cigarstring:
                continue            # skip not clipped read

            m_len = [int(x) for x in re.findall('(\d+)M', read.cigarstring)]    # collect matched lengths
            m_len = max(m_len)                                # leave only longest match from read
            self.match_lengths.append(m_len)                  # add lengths to collection

#           for i in range(read.cigarstring.count('S')):
            boundaries = self._clboundaries(read)
            for cl_seg in boundaries:
                if (cl_seg[2] == "left" and edge == "start") or (cl_seg[2] == "right" and edge == "stop"):
                    clip_temp = self._cltbl_init()
                    clip_temp.at[0, 'ID'] = self._index
                    clip_temp.at[0, 'IS name'] = is_name
                    clip_temp.at[0, 'IS chrom'] = chrom
                    clip_temp.at[0, 'Read name'] = read.query_name
                    clip_temp.at[0, 'left pos'] = cl_seg[0]
                    clip_temp.at[0, 'right pos'] = cl_seg[1]
                    clip_temp.at[0, 'clip position'] = cl_seg[2]
                    clip_temp.at[0, 'junction in read'] = cl_seg[3]
                    if read.is_reverse:
                        clip_temp.at[0, 'reverse'] = True
                    else:
                        clip_temp.at[0, 'reverse'] = False
                    clip_temp.at[0, 'sequence'] = self._clipped_seq(read, cl_seg[0], cl_seg[1])
                    self.clipped_reads = self.clipped_reads.append(clip_temp)
                    self._index = self._index + 1

    # write clipped reads to fasta file with minimum length
    def _write_fasta(self, cl_fasta_name, min_len):
        fasta_file = open(cl_fasta_name, 'w')
        self.clipped_reads.index = self.clipped_reads['ID']
        for index in self.clipped_reads.index:
            if len(self.clipped_reads.at[index, 'sequence']) >= min_len:
                fasta_file.write('>' + str(self.clipped_reads.at[index, 'ID']) + '\n')
                fasta_file.write(str(self.clipped_reads.at[index, 'sequence']) + '\n')
                fasta_file.write('\n')
        fasta_file.close()

    # run blast and write output to xml
    def runblast(self):
        print('Run BLAST for clipped parts of the reads')
        fasta_file = os.path.join(self.outdir, 'cl.fasta')
        blast_out_file = os.path.join(self.outdir, 'cl_blast.out')
        self._write_fasta(fasta_file, self.blast_min)
        blastn_cl = NcbiblastnCommandline(query=fasta_file, db=self.ref_name, evalue=0.001, out=blast_out_file,
                                          outfmt=6, word_size=10)
        blastn_cl()

    # choose left or right coordinate as a clipped junction and orientation relative to junction
    @staticmethod
    def _choosecoord(qleft, qright, lr):
        qcoord = [qleft, qright]
        qorientation = ['left', 'rigth']
        coord = int(qcoord[lr == 'left'])
        orientation = qorientation[not (qcoord[1] > qcoord[0]) ^ (lr == 'left')]
        return coord, orientation

    # parse BALST output
    def parseblast(self):
        print('Collect information from BLAST')
        blast_out = pd.read_csv(os.path.join(self.outdir, 'cl_blast.out'), sep='\t')

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

        blast_out = blast_out[blast_out['pident'] >= 90]  # filter only hits with 90% or higher
        max_in_groups = blast_out.reset_index().groupby(['qseqid'])['bitscore'].max()  # maximum bit score for each read
        count_max = dict()  # for each clipped segment summarize number of best hits
        temp = pd.DataFrame(columns=blast_out.columns)  # temporary dataframe for filtering

        for index in blast_out.index:
            if blast_out.at[index, 'bitscore'] == max_in_groups.at[blast_out.at[index, 'qseqid']]:
                temp = temp.append(blast_out.loc[[index]])
                if blast_out.at[index, 'qseqid'] in count_max.keys():
                    count_max[blast_out.at[index, 'qseqid']] = count_max[blast_out.at[index, 'qseqid']] + 1
                else:
                    count_max[blast_out.at[index, 'qseqid']] = 1

        drop_seg = [cl_seq for cl_seq in count_max.keys() if count_max[cl_seq] > 1]  # collect IDs to drop from the table
        temp = temp[~temp['qseqid'].isin(drop_seg)]  # drop from the table
        for index in temp.index:
            pos, orient = self._choosecoord(temp.at[index, 'sstart'],
                                            temp.at[index, 'send'],
                                            self.clipped_reads.at[temp.at[index, 'qseqid'], 'clip position'])
            temp.at[index, 'pos in ref'] = pos
            temp.at[index, 'orientation'] = orient

        # if temp has any entries
        if temp.size:
            temp['pos in ref'] = temp['pos in ref'].astype(int)
        self.blastout_filtered = self.blastout_filtered.append(temp)

    # create table for junctions description
    def call_junctions(self):
        print('Create junction table')
        for index in self.blastout_filtered.index:
            jtemp = self._jtbl_init()
            read_id = self.blastout_filtered['qseqid'][index]
            jtemp.at[index, 'ID'] = read_id
            jtemp.at[index, 'IS name'] = self.clipped_reads['IS name'][read_id]
            jtemp.at[index, 'IS chrom'] = self.clipped_reads['IS chrom'][read_id]
            jtemp.at[index, 'IS pos'] = self.clipped_reads['clip position'][read_id]
            jtemp.at[index, 'Read name'] = self.clipped_reads['Read name'][read_id]
            jtemp.at[index, 'Chrom'] = self.blastout_filtered['sseqid'][index]
            jtemp.at[index, 'Position'] = self.blastout_filtered['pos in ref'][index]
            jtemp.at[index, 'Orientation'] = self.blastout_filtered['orientation'][index]
            jtemp.at[index, 'Note'] = '-'
            jtemp.at[index, 'Locus tag'] = self.gff.gff_pos[jtemp.at[index, 'Chrom']][jtemp.at[index, 'Position']][0]
            jtemp.at[index, 'Gene'] = self.gff.gff_pos[jtemp.at[index, 'Chrom']][jtemp.at[index, 'Position']][1]
            for b in self.boundaries:
                if b[0] * 2 <= jtemp.at[index, 'Position'] <= b[1] * 2:        # use doubled boundaries
                    jtemp.at[index, 'Note'] = 'IS element'
                    break
            self.junctions = self.junctions.append(jtemp, sort=False)

        self.junctions = self.junctions.reset_index()

    # create summary table that specifies range of junctions in the gene
    def _junction_range(self):
        j_range = self.junctions[['IS name', 'Locus tag', 'Position']].\
            groupby(['IS name', 'Locus tag'], as_index=False).\
            agg({'Position': [min, max]}).\
            rename(columns={'Locus tag': 'Annotation'})

        j_range.columns = ['IS name', 'Annotation', 'min', 'max']
        j_range = j_range.astype({'min': int, 'max': int})
        return j_range

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
            for ann_id, item in self.gff.ann_pos[chrom].items():       #
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
        #aln_depth = self.aln.count_coverage(chrom, start, stop)
        #depth = sum(map(sum, aln_depth))
        #return depth / len(aln_depth[0])  # average depth of the region
        c = pysamstats.load_coverage(self.aln,
                                     chrom=chrom,
                                     start=start,
                                     end=stop + 1,  # pysamstats does not include stop position
                                     truncate=True,
                                     max_depth=100000,
                                     one_based=True,
                                     pad=True)
        return mean(c.reads_all)

    # create report by IS and region
    def report(self):
        print("Create report table")
        self.min_match = min(self.match_lengths)    # find minimum match length
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
        junction_range = self._junction_range()
        self.report_table = pd.merge(self.report_table, junction_range,
                                                    left_on=['IS Name', 'ann'],
                                                    right_on=['IS name', 'Annotation'],
                                                    how='left')
        self.report_table['Depth'] = self.report_table.apply(
            lambda x: self._av_depth(x['chrom'], x['min'], x['max']),
            axis=1
        )

        self.report_table['Frequency'] = self.report_table.apply(
            lambda x: round((x['count'] / 2 * (1 + self.blast_min / self.av_read_len)) / (x['Depth'] * (1 - self.min_match / self.av_read_len)), 4),
            axis=1
        )

        self.report_table = self.report_table[['IS Name', 'ann', 'chrom', 'start', 'stop', 'Frequency', 'Depth']]
        self.report_table.columns = ['IS Name', 'Annotation', 'Chromosome', 'Start', 'Stop', 'Frequency', 'Depth']

# create Circos files
    def create_circos_files(self):
        print('Create CIRCOS files')
        while not os.path.exists(self.data_folder):
            #self.session_id = self._rand_str(5)
            #self.data_folder = "./data" + self.session_id + "/"
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

        text_regions = list()       # to remove duplicates

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
            if self.report_table.iloc[i]['Frequency'] >= self.cutoff:     # draw only lines with cutoff more than specified
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
                                   self.sum_by_region.iloc[i]['stop'],)           # average depth of the region

            # recalculate junction counts to depth
            h_columns = ['chrom', 'start', 'stop']
            h_columns_is = [x for x in self.is_coords.keys()]

            for h in h_columns_is:
                if depth > 0:
                    if self.sum_by_region.iloc[i][h] / depth / 2 >= self.cutoff:
                        histogram.write(' '.join(self.sum_by_region.iloc[i][h_columns].apply(str)) + ' ' +
                                            ','.join(self.sum_by_region.iloc[i][h_columns_is].apply(lambda x: round(((x / depth / 2) * 100), 2)).apply(str)) + '\n')
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

