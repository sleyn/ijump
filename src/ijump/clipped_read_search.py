# Clipped-read search: turn an alignment into a BLAST hit table of candidate
# junction positions. This is the detection pipeline's first algorithmic
# stage -- find reads clipped near IS element boundaries (or, in precise
# mode's backward pass, near estimated insertion sites in the reference),
# BLAST their clipped segments, and keep the best hit per read.
#
# direction: 1 => IS->Ref forward search, 0 => Ref->IS backward search
# (precise mode only). Kept as a single flag rather than split into two
# functions since the two searches share almost all of their logic.
#
# NoInsertionsFound lives here (not in isclipped.py) so this module has no
# import-time dependency on isclipped -- isclipped.py imports it from here
# instead, the same direction it already imports `search`.

import logging
import os
import re
import subprocess
from dataclasses import dataclass, field

import pandas as pd

# Minimum length of a clipped segment to write into the BLAST query FASTA.
# Never overridden by a caller.
BLAST_MIN = 10
# Minimum percent identity of a clipped-read BLAST hit to accept it.
# Never overridden by a caller.
BLAST_MIN_IDENT = 75


class NoInsertionsFound(Exception):
    """The analysis ran correctly and there is nothing to report."""


@dataclass
class Boundary:
    """One area to search for clipped reads.

    For the backward (direction=0) search, ``edge`` and ``is_name`` carry no
    per-element meaning and are always ``"-"`` -- see
    ``ISClipped.set_is_boundaries`` (forward) and the ``backward_boundaries``
    list comprehension in ``ISClipped.run`` (backward) for how each is built.
    """

    start: int
    stop: int
    edge: str
    is_name: str
    chrom: str


@dataclass
class SearchResult:
    clipped_reads: pd.DataFrame
    blast_hits: pd.DataFrame
    match_lengths: list = field(default_factory=list)  # forward (direction=1) only; [] on backward
    read_lengths: int = 0  # forward only; 0 on backward
    n_reads_analyzed: int = 0  # forward only; 0 on backward
    cl_read_cov_overlap: dict = field(
        default_factory=dict
    )  # backward (direction=0) only; {} on forward


# Writes out_file as a side effect. This is the one seam `search` lets tests
# substitute a fake for, so the real blastn subprocess is never invoked by a
# test.
def run_blast_subprocess(query_file, ref_name, out_file):
    logging.info("Run BLAST for clipped parts of the reads")
    try:
        subprocess.run(
            [
                "blastn",
                "-query",
                query_file,
                "-db",
                ref_name,
                "-evalue",
                "0.001",
                "-out",
                out_file,
                "-outfmt",
                "6",
                "-word_size",
                "10",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        logging.error("blastn not found. Is BLAST+ installed and on PATH?")
        raise
    except subprocess.CalledProcessError as e:
        logging.error(f"blastn failed (exit {e.returncode}): {e.stderr.decode(errors='replace')}")
        raise


# For a clipped segment return left, right positions, junction side, coordinate of
# adjacent non-clipped nucleotide.
def _clboundaries(read):
    positions = read.get_reference_positions(full_length=True)
    clipped_segments = list()
    start_clipped_segment = 0
    is_cl_prev = False
    # The segment is "left" if the unmapped part of the read is not aligned at
    # the current position.
    is_left = None
    clipped_part = ""
    junction_pos = 0

    # A read can carry more than one clipped segment (e.g. CIGAR: 30S90M30S).
    for pos_index in range(len(positions)):
        if positions[pos_index] is None and is_cl_prev is False:
            if pos_index == 0:
                is_left = True
            else:
                is_left = False
                clipped_part = "right"
                junction_pos = positions[pos_index - 1]
            start_clipped_segment = pos_index + 1
            is_cl_prev = True
        # If current position is (1a) in aligned part of the read or (1b) it is an end of
        # the read and (2) previous position was in unaligned part of the read then (3)
        # assign this position as the end of the clipped part.
        elif (
            isinstance(positions[pos_index], int) or (pos_index + 1) == len(positions)
        ) and is_cl_prev is True:
            end_clipped_segment = pos_index
            if (pos_index + 1) == len(positions):
                end_clipped_segment = pos_index + 1
            is_cl_prev = False
            if is_left is True:
                junction_pos = positions[pos_index]
                clipped_part = "left"

            clipped_segments.append(
                [start_clipped_segment, end_clipped_segment, clipped_part, junction_pos]
            )
    return clipped_segments


# return clipped portion of a read
def _clipped_seq(read, left, right):
    return read.query_sequence[(left - 1) : right]


# Collect information about coverage that comes from clipped reads outside junction position.
# Mutates cl_read_cov_overlap[chrom] in place rather than returning a value.
def _cl_read_cov_overlap(cl_read_cov_overlap, aln_pairs, chrom):
    if len(aln_pairs) < 3:
        return

    read_pos = [a_pair[0] for a_pair in aln_pairs]
    ref_pos = [a_pair[1] for a_pair in aln_pairs]

    for i in read_pos[1:-1]:
        if i is None:
            continue
        else:
            # Skip a position that is itself a junction.
            if ref_pos[i - 1] is None or ref_pos[i + 1] is None:
                continue
            else:
                cl_read_cov_overlap[chrom][ref_pos[i]] = (
                    cl_read_cov_overlap[chrom].get(ref_pos[i], 0) + 1
                )


# Collect clipped reads from the intervals that do not cross boundaries of a contig.
# The more correct way to collect junction positions would be to find another part of the
# clipped read in the alignment and take its coordinates. CIGAR strings for both parts could
# be not mirrored due to some short repeats (1+nt size) near junction positions.
# direction: 1 => IS->Ref, 0 => Ref->IS (in precise pipeline)
#
# clipped_reads, cl_read_cov_overlap and match_lengths are built locally,
# scoped to this one boundary, and returned rather than mutated on a shared
# object -- `search` (the only caller) merges each call's return into its own
# running totals. `_cl_read_cov_overlap` still mutates the dict it is handed,
# but that dict is this function's own local, never a reference `search`
# holds onto.
def _crtable_ungapped(
    aln,
    index,
    read_lengths,
    n_reads_analyzed,
    boundary,
    direction,
):
    clipped_reads = {}
    cl_read_cov_overlap = {}
    match_lengths = []

    # One is added to convert from 0-based to 1-based system
    for read in aln.fetch(boundary.chrom, boundary.start + 1, boundary.stop + 1):
        if direction:
            if read.infer_read_length():
                read_lengths += read.infer_read_length()
                n_reads_analyzed += 1

        if read.is_unmapped:
            continue

        if "S" not in read.cigarstring:
            continue

        if direction:
            # Longest matched segment, used to calculate the correction
            # coefficient in region_summary/frequency_estimation.
            m_len = [int(x) for x in re.findall(r"(\d+)M", read.cigarstring)]
            m_len = max(m_len)
            match_lengths.append(m_len)
        else:
            # _cl_read_cov_overlap indexes straight into cl_read_cov_overlap[chrom],
            # so that key must exist before the first read seen for a given chrom.
            cl_read_cov_overlap.setdefault(read.reference_name, {})
            _cl_read_cov_overlap(cl_read_cov_overlap, read.aligned_pairs, read.reference_name)

        clipped_segments = _clboundaries(read)
        for cl_seg in clipped_segments:
            # On the IS->Ref search check if read was collected on the correct side of
            # the IS element
            if direction and not (
                (cl_seg[2] == "left" and boundary.edge == "start")
                or (cl_seg[2] == "right" and boundary.edge == "stop")
            ):
                continue

            clip_temp = {
                "ID": index,
                "IS name": boundary.is_name if direction else "-",
                # Contig where IS element is located for IS->Ref search and clipped read of Ref->IS
                "IS_chrom": boundary.chrom,
                "Read name": read.query_name,
                "left pos": cl_seg[0],
                "right pos": cl_seg[1],
                "clip_position": cl_seg[2],
                # Coordinate of junction nucleotide on contig.
                # At IS side for IS->Ref search and for contig for Ref->IS search
                "junction_in_read": cl_seg[3],
                "reverse": True if read.is_reverse else False,
                "sequence": _clipped_seq(read, cl_seg[0], cl_seg[1]),
            }

            # Built as a dict-of-dicts and turned into a DataFrame once at the
            # end (in `search`), since the number of clipped segments isn't
            # known up front and appending row-by-row is far slower.
            clipped_reads[index] = clip_temp

            index = index + 1

    return index, read_lengths, n_reads_analyzed, clipped_reads, cl_read_cov_overlap, match_lengths


# Write clipped parts of reads to FASTA file. Use only parts >= min_len.
def _write_cl_fasta(cl_table, cl_fasta_name, min_len):
    with open(cl_fasta_name, "w") as fasta_file:
        cl_table.index = cl_table["ID"]
        for index in cl_table.index:
            if len(cl_table.at[index, "sequence"]) >= min_len:
                fasta_file.write(">" + str(cl_table.at[index, "ID"]) + "\n")
                fasta_file.write(str(cl_table.at[index, "sequence"]) + "\n")
                fasta_file.write("\n")


# Choose left or right coordinate as a clipped junction and orientation relative to junction.
def _choosecoord(qleft, qright, lr):
    qcoord = [qleft, qright]
    qorientation = ["left", "right"]
    coord = int(qcoord[lr == "left"])
    orientation = qorientation[not (qcoord[1] > qcoord[0]) ^ (lr == "left")]
    return coord, orientation


# Parse BLAST output.
# direction: 1=>IS->Ref, 0=>Ref->IS. For Ref->IS we don't need to remove duplicates,
# we need only one of them
def _parseblast(blast_out_path, direction, cl_table):
    logging.info("Collect information from BLAST")
    if os.path.isfile(blast_out_path) and os.path.getsize(blast_out_path) > 0:
        blast_out = pd.read_csv(blast_out_path, sep="\t")
    else:
        logging.info("No BLAST hits were found.")
        raise NoInsertionsFound("No BLAST hits were found.")

    blast_out.columns = [
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
    ]

    # Filter only hits with identity [BLAST_MIN_IDENT]% or higher. Default: 75%.
    blast_out = blast_out[blast_out["pident"] >= BLAST_MIN_IDENT]

    idx_max = blast_out.groupby("qseqid")["bitscore"].transform(max) == blast_out["bitscore"]
    temp = blast_out[idx_max].copy()

    if direction:
        # A read with more than one equally-best hit is ambiguous and dropped.
        temp["count"] = temp.groupby("qseqid")["qseqid"].transform("count")
        temp = temp[temp["count"] == 1]
        temp = temp.drop(columns=["count"])
    else:
        # Ties are fine here -- just pick one deterministically.
        temp["rank"] = temp.groupby("qseqid")["bitscore"].rank(method="first", ascending=True)
        temp = temp[temp["rank"] == 1]
        temp = temp.drop(columns=["rank"])

    for index in temp.index:
        pos, orient = _choosecoord(
            temp.at[index, "sstart"],
            temp.at[index, "send"],
            cl_table.at[temp.at[index, "qseqid"], "clip_position"],
        )
        temp.at[index, "pos_in_ref"] = pos
        temp.at[index, "orientation"] = orient

    if temp.size:
        temp["pos_in_ref"] = temp["pos_in_ref"].astype(int)
    else:
        logging.info("No significant BLAST hits.")
        raise NoInsertionsFound("No significant BLAST hits.")

    return temp


# Turn an alignment into a BLAST hit table of candidate junction positions.
#
# `boundaries` is a list of `Boundary` entries -- the shape
# ISClipped.set_is_boundaries builds for the forward search. For the backward
# search, the caller builds the same shape from the reference-regions table,
# with edge=is_name='-'.
#
# Raises NoInsertionsFound (without ever calling `run_blast`) if no clipped
# reads are collected.
def search(aln, boundaries, ref_name, workdir, direction, run_blast=run_blast_subprocess):
    clipped_reads_dict = {}
    match_lengths = []
    read_lengths = 0
    n_reads_analyzed = 0
    cl_read_cov_overlap = {ref: {} for ref in aln.references}
    index = 0

    for boundary in boundaries:
        logging.info(
            "Collect clipped reads: "
            + " ".join(
                str(x)
                for x in (
                    boundary.start,
                    boundary.stop,
                    boundary.edge,
                    boundary.is_name,
                    boundary.chrom,
                )
            )
        )

        (
            index,
            read_lengths,
            n_reads_analyzed,
            boundary_clipped_reads,
            boundary_cl_read_cov_overlap,
            boundary_match_lengths,
        ) = _crtable_ungapped(aln, index, read_lengths, n_reads_analyzed, boundary, direction)

        # Merge this boundary's contribution into the running totals. Index
        # ranges never overlap between boundaries (each call picks up where
        # the previous left off), so clipped_reads_dict can just be updated;
        # cl_read_cov_overlap positions can recur across boundaries and must
        # be summed, not overwritten.
        clipped_reads_dict.update(boundary_clipped_reads)
        match_lengths.extend(boundary_match_lengths)
        for chrom, positions in boundary_cl_read_cov_overlap.items():
            chrom_cov = cl_read_cov_overlap.setdefault(chrom, {})
            for pos, count in positions.items():
                chrom_cov[pos] = chrom_cov.get(pos, 0) + count

    clipped_reads = pd.DataFrame.from_dict(clipped_reads_dict, "index")

    if clipped_reads.size == 0:
        message = (
            "No clipped reads were found."
            if direction
            else "No clipped reads were found near estimated insertion sites."
        )
        raise NoInsertionsFound(message)

    if direction:
        fasta_file = os.path.join(workdir, "cl.fasta")
        blast_out_file = os.path.join(workdir, "cl_blast.out")
    else:
        fasta_file = os.path.join(workdir, "cl_bwrd.fasta")
        blast_out_file = os.path.join(workdir, "cl_blast_bwds.out")

    _write_cl_fasta(clipped_reads, fasta_file, BLAST_MIN)
    run_blast(fasta_file, ref_name, blast_out_file)
    blast_hits = _parseblast(blast_out_file, direction, clipped_reads)

    if direction:
        return SearchResult(
            clipped_reads=clipped_reads,
            blast_hits=blast_hits,
            match_lengths=match_lengths,
            read_lengths=read_lengths,
            n_reads_analyzed=n_reads_analyzed,
            cl_read_cov_overlap={},
        )
    else:
        return SearchResult(
            clipped_reads=clipped_reads,
            blast_hits=blast_hits,
            match_lengths=[],
            read_lengths=0,
            n_reads_analyzed=0,
            cl_read_cov_overlap=cl_read_cov_overlap,
        )
