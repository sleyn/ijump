#!/usr/bin/env python3

import argparse
from os.path import join as join_path

import pandas as pd

from . import is_annotation, is_table


def main():
    parser = argparse.ArgumentParser(description="Parse BLAST output Genome vs ISFinder.")
    parser.add_argument(
        "-b",
        "--blast_out",
        type=str,
        action="store",
        help="BLAST output for parsing. Require outfmt 6.",
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        action="store",
        required=True,
        help="Reference FASTA the BLAST search was run against. Each called element "
        "is extracted from it and compared with all the others, so copies of one "
        "mobile element share a cluster.",
    )
    is_annotation.add_cluster_arguments(parser)
    # parser.add_argument('-c', '--csv', type=str, action='store',
    #                      help='CSV description of IS Finder database mobile elements.')
    parser.add_argument("-o", "--outdir", type=str, default=".", help="Output directory")
    args = parser.parse_args()

    # Outfmt 6 BLAST output
    blast_out = pd.read_csv(
        args.blast_out,
        sep="\t",
        names=[
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
        ],
    )

    # Filter high e-value
    blast_out = blast_out[blast_out.evalue <= 1e-30]
    blast_out = blast_out.sort_values(by=["bitscore"], ascending=False)

    # A subject id is name_family_group -- keep all three, splitting from the
    # right so names carrying an underscore of their own survive whole.
    blast_out[["sseqid", "family", "group"]] = pd.DataFrame(
        blast_out["sseqid"].apply(is_table.parse_subject_id).tolist(),
        columns=["sseqid", "family", "group"],
        index=blast_out.index,
    )

    # Information about ISFinder databasse records
    # db_info = pd.read_csv(args.csv, usecols=['Name', 'Family', 'Group', 'Origin', 'Length'])

    # check if IS elements are overlapping
    occupied_pos = {}

    contigs_list = blast_out.qseqid.unique()
    for contig in contigs_list:
        occupied_pos[contig] = set()

    # check_overlap: True - no overlap, False - overlap
    check_overlap = []

    for is_element in blast_out.itertuples(name=None, index=False):
        contig = is_element[0]
        qstart = is_element[6]
        qend = is_element[7]
        overlap = occupied_pos[contig].intersection(range(qstart, qend + 1))

        # check if IS element occupies a new region at least 75% of its length
        if len(overlap) / (qend - qstart + 1) >= 0.75:
            # Set check to remove
            check_overlap.append(False)
        else:
            # Set check to keep
            check_overlap.append(True)
            # add positions occupied by IS element
            occupied_pos[contig].update(range(qstart, qend + 1))

    # drop IS elements that did not pass overlap check
    blast_out = blast_out[check_overlap]
    blast_out["Row_num"] = blast_out.groupby(["sseqid"]).cumcount() + 1
    blast_out["sseqid"] = blast_out.apply(lambda x: x.sseqid + "_" + str(x.Row_num), axis=1)

    blast_out = blast_out.rename(
        columns={
            "sseqid": "is_name",
            "qseqid": "contig",
            "qstart": "start",
            "qend": "stop",
        }
    )

    blast_out = is_annotation.annotate_and_cluster(
        blast_out,
        args.ref,
        identity=args.cluster_identity,
        coverage=args.cluster_coverage,
    )

    is_table.write_is_table(blast_out, join_path(args.outdir, "ISTable_processing.txt"))


if __name__ == "__main__":
    main()
