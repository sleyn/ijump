# Clipped-read search: turn an alignment into a BLAST hit table of candidate
# junction positions. This is the detection pipeline's first algorithmic
# stage -- find reads clipped near IS element boundaries (or, in precise
# mode's backward pass, near estimated insertion sites in the reference),
# BLAST their clipped segments, and keep the best hit per read.
#
# direction: 1 => IS->Ref forward search, 0 => Ref->IS backward search
# (precise mode only). Kept as a single flag rather than split into two
# functions -- see ticket 10 ("Out of scope").
#
# NoInsertionsFound lives here (not in isclipped.py) so this module has no
# import-time dependency on isclipped -- isclipped.py imports it from here
# instead, the same direction it already imports `search`.

import os
import re
import logging
from dataclasses import dataclass, field

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline

# Minimum length of a clipped segment to write into the BLAST query FASTA.
# Was ISClipped.blast_min; never overridden by a caller.
BLAST_MIN = 10
# Minimum percent identity of a clipped-read BLAST hit to accept it.
# Was ISClipped.blast_min_ident; never overridden by a caller.
BLAST_MIN_IDENT = 75


class NoInsertionsFound(Exception):
    """The analysis ran correctly and there is nothing to report."""


@dataclass
class SearchResult:
    clipped_reads: pd.DataFrame        # replaces self.clipped_reads / self.clipped_reads_bwrd
    blast_hits: pd.DataFrame           # replaces self.blastout_filtered
    match_lengths: list = field(default_factory=list)  # forward (direction=1) only; [] on backward
    read_lengths: int = 0              # forward only; 0 on backward
    n_reads_analyzed: int = 0          # forward only; 0 on backward
    cl_read_cov_overlap: dict = field(default_factory=dict)  # backward (direction=0) only; {} on forward


# Run blastn on a FASTA of clipped-read segments against the reference.
# Writes out_file as a side effect, same as today's ISClipped.runblast
# -- the one seam `search` lets tests substitute a
# fake for, so the real blastn subprocess is never invoked by a test.
def run_blast_subprocess(query_file, ref_name, out_file):
    logging.info('Run BLAST for clipped parts of the reads')
    blastn_cl = NcbiblastnCommandline(query=query_file, db=ref_name, evalue=0.001, out=out_file,
                                      outfmt=6, word_size=10)
    blastn_cl()


# For a clipped segment return left, right positions, junction side, coordinate of adjacent non-clipped nucleotide.
# Moved verbatim from ISClipped._clboundaries.
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
# Moved verbatim from ISClipped._clipped_seq. Not one
# of the eight functions ticket 10 names explicitly, but _crtable_ungapped
# depends on it, so it moves along with its caller.
def _clipped_seq(read, left, right):
    return read.query_sequence[(left - 1):right]


# Collect information about coverage that comes from clipped reads outside junction position.
# Moved from ISClipped._cl_read_cov_overlap: the
# `self.cl_read_cov_overlap[chrom]` mutation becomes a `cl_read_cov_overlap`
# parameter mutated in place (dicts are mutable, so no return value is
# needed -- same pattern as the original, which also discarded its `return 0`
# early-exit's value).
def _cl_read_cov_overlap(cl_read_cov_overlap, aln_pairs, chrom):
    if len(aln_pairs) < 3:
        return

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
                cl_read_cov_overlap[chrom][ref_pos[i]] = cl_read_cov_overlap[chrom].get(ref_pos[i], 0) + 1


# Collect clipped reads from the intervals that do not cross boundaries of a contig.
# The more correct way to collect junction positions would be to find another part of the
# clipped read in the alignment and take its coordinates. CIGAR strings for both parts could
# be not mirrored due to some short repeats (1+nt size) near junction positions.
# direction: 1 => IS->Ref, 0 => Ref->IS (in precise pipeline)
#
# Moved from ISClipped._crtable_ungapped. `self._index`
# becomes the `index` parameter/return value; `self.match_lengths` and
# `self.cl_read_cov_overlap`/`clipped_reads_dict` become parameters mutated
# in place; `self.read_lengths`/`self.n_reads_analyzed` become
# parameters/return values (plain ints can't be mutated through a reference).
def _crtable_ungapped(aln, clipped_reads_dict, cl_read_cov_overlap, match_lengths,
                       index, read_lengths, n_reads_analyzed,
                       chrom, start, stop, edge, is_name, direction):  # generate clipped read table
    # One is added to convert from 0-based to 1-based system
    for read in aln.fetch(chrom, start + 1, stop + 1):
        # Add read length to collection of lengths.
        if direction:
            if read.infer_read_length():
                read_lengths += read.infer_read_length()
                n_reads_analyzed += 1

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
            match_lengths.append(m_len)
        else:
            # If it is Ref->IS direction of search:
            # Add coverage from aligned positions of clipped reads that are not junctions.
            _cl_read_cov_overlap(cl_read_cov_overlap, read.aligned_pairs, read.reference_name)

        # Get clipped segments coordinates from the read
        boundaries = _clboundaries(read)
        for cl_seg in boundaries:
            # On the IS->Ref search check if read was collected on the correct side of the IS element
            if direction and \
                    not ((cl_seg[2] == "left" and edge == "start") or (cl_seg[2] == "right" and edge == "stop")):
                continue

            clip_temp = {
                # Unique read ID
                'ID': index,
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
                'sequence': _clipped_seq(read, cl_seg[0], cl_seg[1])
            }

            # Add clipped read information to dictionary to build DataFrame.
            # It is faster than append segments-by-segments to the existing DataFrame.
            # As we don't know number of clipped segments we could not create and empty DataFrame of required size.
            clipped_reads_dict[index] = clip_temp

            index = index + 1

    return index, read_lengths, n_reads_analyzed


# Write clipped parts of reads to FASTA file. Use only parts >= min_len.
# Moved from ISClipped._write_cl_fasta. The
# direction-based table selection moves to the call site in `search`, since
# each `search` call already only ever has one table in scope.
def _write_cl_fasta(cl_table, cl_fasta_name, min_len):
    with open(cl_fasta_name, 'w') as fasta_file:
        cl_table.index = cl_table['ID']
        for index in cl_table.index:
            if len(cl_table.at[index, 'sequence']) >= min_len:
                fasta_file.write('>' + str(cl_table.at[index, 'ID']) + '\n')
                fasta_file.write(str(cl_table.at[index, 'sequence']) + '\n')
                fasta_file.write('\n')


# Choose left or right coordinate as a clipped junction and orientation relative to junction.
# Moved verbatim from ISClipped._choosecoord.
def _choosecoord(qleft, qright, lr):
    qcoord = [qleft, qright]
    qorientation = ['left', 'right']
    coord = int(qcoord[lr == 'left'])
    orientation = qorientation[not (qcoord[1] > qcoord[0]) ^ (lr == 'left')]
    return coord, orientation


# Parse BLAST output.
# direction: 1=>IS->Ref, 0=>Ref->IS. For Ref->IS we don't need to remove duplicates, we need only one of them
#
# Moved from ISClipped.parseblast. `self.blastout_filtered`
# assignment becomes a return value; `self.clipped_reads`/`self.clipped_reads_bwrd`
# become the `cl_table` parameter (chosen at the call site, same as
# _write_cl_fasta above); `self.blast_min_ident` becomes the BLAST_MIN_IDENT
# module constant (see its definition above).
def _parseblast(blast_out_path, direction, cl_table):
    logging.info('Collect information from BLAST')
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

    # Filter only hits with identity [BLAST_MIN_IDENT]% or higher. Default: 75%.
    blast_out = blast_out[blast_out['pident'] >= BLAST_MIN_IDENT]

    idx_max = blast_out.groupby('qseqid')['bitscore'].transform(max) == blast_out['bitscore']
    # Temporary dataframe for filtering with only best hits by bitscore.
    temp = blast_out[idx_max].copy()

    if direction:
        temp['count'] = temp.groupby('qseqid')['qseqid'].transform('count')
        # Leave only hits with one best hit.
        temp = temp[temp['count'] == 1]
        temp = temp.drop(columns=['count'])
    else:
        temp['rank'] = temp.groupby('qseqid')['bitscore'].rank(method="first", ascending=True)
        # Leave only first best hit.
        temp = temp[temp['rank'] == 1]
        temp = temp.drop(columns=['rank'])

    for index in temp.index:
        pos, orient = _choosecoord(temp.at[index, 'sstart'],
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

    return temp


# Turn an alignment into a BLAST hit table of candidate junction positions.
#
# `boundaries` is a list of [start, stop, edge, is_name, chrom] entries --
# the shape ISClipped.boundaries already has for the forward search
# (ISClipped.set_is_boundaries). For the backward search, the caller
# builds the same shape from the reference-regions table (what
# ISClipped.crtable_bwds used to do internally with edge=is_name='-').
#
# Raises NoInsertionsFound (without ever calling `run_blast`) if no clipped
# reads are collected -- same short-circuit today's ISClipped.run() gets
# from calling check_data_presence_in_df between collect_clipped_reads()/
# crtable_bwds() and runblast().
def search(aln, boundaries, ref_name, workdir, direction, run_blast=run_blast_subprocess):
    clipped_reads_dict = {}
    match_lengths = []
    read_lengths = 0
    n_reads_analyzed = 0
    cl_read_cov_overlap = {ref: {} for ref in aln.references}
    index = 0

    for b in boundaries:
        logging.info('Collect clipped reads: ' + ' '.join(str(x) for x in b))

        start_collection, stop_collection, edge_of_is, name_of_is, chrom = b[0], b[1], b[2], b[3], b[4]

        index, read_lengths, n_reads_analyzed = _crtable_ungapped(
            aln, clipped_reads_dict, cl_read_cov_overlap, match_lengths,
            index, read_lengths, n_reads_analyzed,
            chrom, start_collection, stop_collection, edge_of_is, name_of_is, direction,
        )

    clipped_reads = pd.DataFrame.from_dict(clipped_reads_dict, "index")

    if clipped_reads.size == 0:
        message = 'No clipped reads were found.' if direction \
            else 'No clipped reads were found near estimated insertion sites.'
        raise NoInsertionsFound(message)

    if direction:
        fasta_file = os.path.join(workdir, 'cl.fasta')
        blast_out_file = os.path.join(workdir, 'cl_blast.out')
    else:
        fasta_file = os.path.join(workdir, 'cl_bwrd.fasta')
        blast_out_file = os.path.join(workdir, 'cl_blast_bwds.out')

    _write_cl_fasta(clipped_reads, fasta_file, BLAST_MIN)
    run_blast(fasta_file, ref_name, blast_out_file)
    blast_hits = _parseblast(blast_out_file, direction, clipped_reads)

    if direction:
        return SearchResult(
            clipped_reads=clipped_reads, blast_hits=blast_hits,
            match_lengths=match_lengths, read_lengths=read_lengths,
            n_reads_analyzed=n_reads_analyzed, cl_read_cov_overlap={},
        )
    else:
        return SearchResult(
            clipped_reads=clipped_reads, blast_hits=blast_hits,
            match_lengths=[], read_lengths=0,
            n_reads_analyzed=0, cl_read_cov_overlap=cl_read_cov_overlap,
        )
