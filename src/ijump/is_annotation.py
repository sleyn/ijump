"""The annotate-and-cluster core every IS-table back-end runs.

A back-end's job is to produce the four locus columns -- name, contig, start,
stop -- from whatever it reads: an ISFinder BLAST search of a genome, an
operator's hand-curated table, an ISEScan run. What happens *after* that is the
same work in every case, and it is work with rules in it: which loci are copies
of one mobile element (``is_clustering``), and which pairs are one copy the
assembler cut at a contig seam (``origin_spanning``). Those rules live here so
the back-ends share them rather than each carrying its own copy to drift.
"""

import os
from typing import Union

import pandas as pd

from ijump import is_clustering, origin_spanning

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]


def annotate_and_cluster(
    table: pd.DataFrame,
    reference: Path,
    identity: float = is_clustering.IDENTITY_DEFAULT,
    coverage: float = is_clustering.COVERAGE_DEFAULT,
) -> pd.DataFrame:
    """Fill ``cluster``, ``wraps_origin`` and ``element_id`` on a locus table.

    ``table`` carries at least ``is_name``, ``contig``, ``start`` and ``stop``;
    it is not modified, and the returned frame is a copy.
    """
    annotated = table.copy()

    # Which rows are copies of one mobile element is a question about their
    # sequences, not about their names: a name is the nearest database entry, and
    # two fragments of one element land on different entries while two distinct
    # elements can land on the same one.
    annotated["cluster"] = is_clustering.annotate(
        annotated, reference, identity=identity, coverage=coverage
    )

    # A copy the assembler cut in half at a contig seam is two rows of the table.
    # Clustering has already united them; this says so on the rows themselves.
    annotated[["wraps_origin", "element_id"]] = origin_spanning.annotate(annotated, reference)

    return annotated
