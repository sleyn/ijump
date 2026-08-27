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
        help="Reference FASTA the BLAST search was run against. Each called locus "
        "is extracted from it and compared with all the others, so copies of one "
        "mobile element share a cluster.",
    )
    is_annotation.add_cluster_arguments(parser)
    parser.add_argument("-o", "--outdir", type=str, default=".", help="Output directory")
    args = parser.parse_args()

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

    blast_out = blast_out[blast_out.evalue <= 1e-30]
    blast_out = blast_out.sort_values(by=["bitscore"], ascending=False)

    # A subject id is name_family_group -- keep all three, splitting from the
    # right so names carrying an underscore of their own survive whole.
    blast_out[["sseqid", "family", "group"]] = pd.DataFrame(
        blast_out["sseqid"].apply(is_table.parse_subject_id).tolist(),
        columns=["sseqid", "family", "group"],
        index=blast_out.index,
    )

    occupied_pos = {}

    contigs_list = blast_out.qseqid.unique()
    for contig in contigs_list:
        occupied_pos[contig] = set()

    # True keeps a hit, False drops it as an overlap.
    check_overlap = []

    for is_element in blast_out.itertuples(name=None, index=False):
        contig = is_element[0]
        qstart = is_element[6]
        qend = is_element[7]
        overlap = occupied_pos[contig].intersection(range(qstart, qend + 1))

        if len(overlap) / (qend - qstart + 1) >= 0.75:
            check_overlap.append(False)
        else:
            check_overlap.append(True)
            occupied_pos[contig].update(range(qstart, qend + 1))

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
