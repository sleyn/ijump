#!/usr/bin/env python3
"""Turn ISEScan's results into an IS table.

iJump reads ISEScan's output; it never runs ISEScan. That was a cost decision:
ISEScan is x64-only (needing emulation on arm64), wants a library symlink
workaround, and takes about fourteen minutes on the reference genome -- in a
project whose runtime dependencies already cannot be resolved by the standard
installer. Reading its output puts that burden only on operators who choose it.

Why support it at all: the two annotations are complementary. On the reference
genome ISEScan finds three copies of an element with no ISFinder database hit
whatsoever, so an element absent from ISFinder and actively jumping is invisible
to iJump otherwise. In the other direction ISEScan needs terminal repeats plus an
ORF, so the 76 bp ``IS17`` remnant is structurally invisible to it. Running both
and taking the union would be the best annotation available for this genome, but
the two disagree on span where both fire -- one locus is 977 bp to one and
2299 bp to the other, and span drives the boundary search windows -- so the union
needs a documented rule for that conflict. That is a scientific judgement and is
deliberately not made here.

Everything after the locus columns is ``is_annotation.annotate_and_cluster``,
shared with the other back-ends.
"""

import argparse
import os
from os.path import join as join_path
from typing import Union

import pandas as pd

from ijump import is_annotation, is_clustering, is_table

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# The ISEScan `.tsv` columns this reads. ISEScan writes many more; these are the
# ones an IS table needs.
SEQ_ID = "seqID"
FAMILY = "family"
CLUSTER = "cluster"
IS_BEGIN = "isBegin"
IS_END = "isEnd"

REQUIRED_COLUMNS = (SEQ_ID, FAMILY, CLUSTER, IS_BEGIN, IS_END)


class NotISEScanOutput(Exception):
    """A file handed in as ISEScan results does not have its columns."""


def read_isescan(path: Path) -> pd.DataFrame:
    """Read ISEScan's tab-separated results as the four locus columns.

    Every call becomes a row. This back-end converts; deciding which calls to
    believe is the operator's, made before the file reaches here.
    """
    results = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    missing = [column for column in REQUIRED_COLUMNS if column not in results.columns]
    if missing:
        raise NotISEScanOutput(
            f"{os.fspath(path)} is missing the ISEScan column(s): {', '.join(missing)}. "
            "This reads ISEScan's tab-separated `.tsv` results -- not its `.out`, "
            "which is a fixed-width report, nor its `.csv`, which is the same "
            "content comma-separated."
        )

    loci = pd.DataFrame(
        {
            # ISEScan reports no element name, only a family and its own cluster
            # id. The id plus a copy number is the unique per-locus name the IS
            # table wants, in the `<element>_<n>` shape the other back-ends use --
            # so base_is_name strips the counter and leaves the element behind.
            "is_name": _copy_numbered(results[CLUSTER]),
            "contig": results[SEQ_ID],
            "start": results[IS_BEGIN],
            "stop": results[IS_END],
            "family": results[FAMILY],
        }
    )
    # `group` and `pident` stay empty on purpose, not for want of filling in:
    # ISEScan reports its own numeric cluster rather than an ISFinder group, and
    # no identity against any database. Filing the cluster id under `group` would
    # put a different kind of value under a name that already means something.
    # `cluster`, `wraps_origin` and `element_id` are filled by annotate_and_cluster.
    return is_table.with_all_columns(loci).reset_index(drop=True)


def convert(
    path: Path,
    reference: Path,
    identity: float = is_clustering.IDENTITY_DEFAULT,
    coverage: float = is_clustering.COVERAGE_DEFAULT,
) -> pd.DataFrame:
    """Read ISEScan results and annotate them into a finished IS table."""
    return is_annotation.annotate_and_cluster(
        read_isescan(path), reference, identity=identity, coverage=coverage
    )


def _copy_numbered(cluster_ids: pd.Series) -> list:
    """``IS701_225`` three times becomes ``IS701_225_1/_2/_3``, in file order."""
    counts: dict = {}
    names = []
    for cluster_id in cluster_ids:
        counts[cluster_id] = counts.get(cluster_id, 0) + 1
        names.append(f"{cluster_id}_{counts[cluster_id]}")
    return names


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Convert ISEScan results into an iJump IS table. iJump reads "
        "ISEScan's output and never runs ISEScan -- run it yourself and pass the "
        "tab-separated results here.",
    )
    parser.add_argument(
        "-i",
        "--isescan",
        type=str,
        required=True,
        help="ISEScan's tab-separated results (`<genome>.fna.tsv`).",
    )
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        required=True,
        help="Reference FASTA ISEScan was run on. Each locus is extracted from it "
        "so copies of one mobile element share a cluster.",
    )
    is_annotation.add_cluster_arguments(parser)
    parser.add_argument("-o", "--outdir", type=str, default=".", help="Output directory")
    return parser


def main():
    args = build_arg_parser().parse_args()
    converted = convert(
        args.isescan,
        args.ref,
        identity=args.cluster_identity,
        coverage=args.cluster_coverage,
    )
    is_table.write_is_table(converted, join_path(args.outdir, "ISTable_processing.txt"))


if __name__ == "__main__":
    main()
