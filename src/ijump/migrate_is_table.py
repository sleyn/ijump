#!/usr/bin/env python3
"""Bring a legacy four-column IS table up to the current format.

An operator who has curated an IS table by hand should not have to throw it away
because iJump started needing annotation it never carried. This back-end takes
that table, keeps every coordinate exactly as written, and fills in what the
current format wants: family and group, the percent identity behind them, the
cluster, and the origin-spanning flags.

Family and group cannot be recovered from a four-column file -- it never held
them -- so they are re-derived by searching each locus against the ISFinder
database, the same database the primary back-end's names come from. Everything
after that is ``is_annotation.annotate_and_cluster``, shared with the primary
back-end so the two cannot drift apart on what a cluster is.

This is the remedy that ``ijump run`` names when it meets a table with no
cluster column.
"""

import argparse
import logging
import os
import tempfile
from os.path import join as join_path
from typing import Dict, Tuple, Union

import pandas as pd

from ijump import is_annotation, is_clustering, is_table

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# Family, group and percent identity of one locus's best database hit.
Annotation = Tuple[str, str, str]

# A locus that matched nothing in the database. Kept in the table with its
# coordinates intact -- the operator put it there on purpose, and this back-end
# annotates, it does not re-call loci.
UNANNOTATED: Annotation = ("", "", "")


def search_database(sequences: Dict[str, str], database: Path) -> Dict[str, Annotation]:
    """Search each locus against the ISFinder database, keeping its best hit.

    Best by bitscore, which is how the primary back-end picks too.

    The e-value threshold is SENSITIVE_BLASTN's 1e-5, not the 1e-30 the primary
    back-end filters on. That back-end is *calling* loci from a whole genome,
    where a lax threshold invents them; this one is annotating loci an operator
    already decided on, and the shortest of them is 76 bp -- short enough that
    1e-30 is unreachable however good the match. Refusing to annotate a locus
    because it is short would be the wrong answer to a question already settled.

    A locus with no hit comes back UNANNOTATED rather than being dropped.
    """
    if not sequences:
        return {}

    with tempfile.TemporaryDirectory() as workdir:
        query = os.path.join(workdir, "loci.fna")
        with open(query, "w") as fasta:
            for name, sequence in sequences.items():
                fasta.write(f">{name}\n{sequence}\n")

        out_file = os.path.join(workdir, "loci_vs_isfinder.out")
        is_clustering.run_blast_command(
            [
                "blastn",
                # The sensitive task, not the binary's megablast default: the
                # shortest loci are under 100 bp and megablast's seed would miss
                # them entirely. Same reasoning as is_clustering.all_vs_all_search.
                "-task",
                "blastn",
                "-query",
                query,
                "-db",
                os.fspath(database),
                "-out",
                out_file,
                "-outfmt",
                "6 qseqid sseqid pident bitscore",
                "-word_size",
                "11",
                "-dust",
                "no",
                "-evalue",
                "1e-5",
            ]
        )
        return _best_hits(out_file, sequences)


def migrate(
    table: pd.DataFrame,
    reference: Path,
    database: Path,
    identity: float = is_clustering.IDENTITY_DEFAULT,
    coverage: float = is_clustering.COVERAGE_DEFAULT,
) -> pd.DataFrame:
    """Annotate ``table``'s loci, leaving their coordinates untouched."""
    loci = is_clustering.loci_from_table(table)
    sequences = is_clustering.extract_loci(reference, loci)
    annotations = search_database(sequences, database)

    migrated = table.copy()
    for column, position in (("family", 0), ("group", 1), ("pident", 2)):
        migrated[column] = [
            annotations.get(str(name), UNANNOTATED)[position] for name in migrated["is_name"]
        ]

    unannotated = sorted(name for name, hit in annotations.items() if hit == UNANNOTATED)
    if unannotated:
        logging.warning(
            "No database hit for: %s. These loci keep their coordinates and are "
            "clustered with the rest, but carry no family or group.",
            ", ".join(unannotated),
        )

    return is_annotation.annotate_and_cluster(
        migrated, reference, identity=identity, coverage=coverage
    )


def _best_hits(out_file: Path, sequences: Dict[str, str]) -> Dict[str, Annotation]:
    """Pick each locus's highest-bitscore hit and split its subject id."""
    best: Dict[str, Tuple[float, Annotation]] = {}
    with open(out_file) as blast_out:
        for line in blast_out:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            query, subject, pident, bitscore = fields[:4]
            score = float(bitscore)
            if query in best and best[query][0] >= score:
                continue
            _, family, group = is_table.parse_subject_id(subject)
            # Through float, not BLAST's raw text: the primary back-end reads its
            # hits with pd.read_csv, so its pident reaches the table as a float
            # and is written "100.0". Keeping "100.000" here would put the same
            # identity in the table two ways depending on which back-end wrote it.
            best[query] = (score, (family, group, str(float(pident))))

    return {name: best[name][1] if name in best else UNANNOTATED for name in sequences}


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Migrate a legacy four-column IS table to the current format, "
        "preserving its coordinates and adding family, group and cluster annotation.",
    )
    parser.add_argument(
        "-i",
        "--is-table",
        type=str,
        required=True,
        help="IS table to migrate. The legacy four-column format (is_name, contig, "
        "start, stop) or any table the current reader accepts.",
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        required=True,
        help="Reference FASTA the table's coordinates refer to. Each locus is "
        "extracted from it, both to search the database and to cluster.",
    )
    parser.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="ISFinder BLAST database (the -out prefix given to makeblastdb), "
        "searched to recover each locus's family and group.",
    )
    is_annotation.add_cluster_arguments(parser)
    parser.add_argument("-o", "--outdir", type=str, default=".", help="Output directory")
    return parser


def main():
    args = build_arg_parser().parse_args()

    table = is_table.read_is_table(args.is_table)
    migrated = migrate(
        table,
        args.ref,
        args.database,
        identity=args.cluster_identity,
        coverage=args.cluster_coverage,
    )
    is_table.write_is_table(migrated, join_path(args.outdir, "ISTable_processing.txt"))


if __name__ == "__main__":
    main()
