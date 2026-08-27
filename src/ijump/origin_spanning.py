"""Find the IS copies an assembler cut in half at a contig boundary.

A circular replicon has to be broken somewhere to be written as a linear contig,
and nothing stops the break landing inside an IS element. When it does, one copy
is called as two loci: one ending at the last base of the contig, one starting at
its first. ``Test/A_baumannii_assembly.fna`` has exactly that -- ``IS17_1`` covers
the last 144 bases of ``NODE_2`` and ``IS17_2`` its first 76, one ISAba12-like
element with the assembly's seam through the middle.

Both rows are kept, and kept as they are. They are genuinely separate spans, the
boundary search needs both, and joining them into a single ``start > stop`` row
would break every consumer that assumes a row's start precedes its stop --
``set_is_boundaries``, ``region_summary``'s overlap logic -- while claiming a
220 bp element where the truth is a 1039 bp one with a hole in the middle.

So this changes no coordinates and merges nothing. ``is_clustering`` already puts
the two halves in one cluster by transitivity; what is missing is any sign in the
output that the assembler broke a contig inside an IS copy. This module supplies
it: a ``wraps_origin`` flag and an ``element_id`` shared by the halves of one
broken copy.
"""

import os
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple, Union

import pandas as pd
import pysam

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# How close to a boundary a span has to come to count as touching it. The seam
# is not always clean: the reference genome's own case has ``IS17_2`` starting at
# base 2 rather than 1, because the alignment that called it frayed a base short.
# Wide enough for that fraying, far narrower than any IS element, so a locus
# sitting inside the contig cannot reach a boundary by accident.
BOUNDARY_MARGIN = 20


@dataclass(frozen=True)
class Span:
    """One row of the IS table, as boundary detection sees it."""

    name: str
    contig: str
    start: int
    stop: int
    cluster: str


def origin_spanning_elements(
    spans: Sequence[Span],
    contig_lengths: Dict[str, int],
) -> Dict[str, str]:
    """``{locus name: element id}`` for the loci that span a contig origin.

    A locus spans the origin when it shares a contig *and* a cluster with a locus
    at the opposite boundary of that contig. Both halves are required: an IS copy
    at the end of a contig is an everyday event, and only a counterpart at the
    origin of the same contig makes it half of something. Loci that do not span
    the origin are absent from the result rather than mapped to an empty id, so
    the caller can tell the two apart without comparing strings.

    Pure -- the contig lengths are passed in -- so the rule is testable without a
    reference FASTA.
    """
    assigned: Dict[str, str] = {}
    counts: Dict[str, int] = {}

    grouped = _by_contig_and_cluster(spans)
    # Sorted, so the numbering an element id carries depends on the table's
    # content rather than on the order its rows happened to arrive in.
    for contig, cluster in sorted(grouped):
        members = grouped[(contig, cluster)]
        length = contig_lengths.get(contig)
        if length is None:
            raise ValueError(
                f"contig {contig!r} of IS element {members[0].name!r} has no length -- "
                "is this the reference the IS table was called against?"
            )

        at_origin = [span for span in members if span.start <= 1 + BOUNDARY_MARGIN]
        at_end = [span for span in members if span.stop >= length - BOUNDARY_MARGIN]
        halves = {span.name: span for span in at_origin + at_end}
        # One locus touching both boundaries is a contig shorter than the element,
        # not an element broken across the seam.
        if not at_origin or not at_end or len(halves) < 2:
            continue

        counts[cluster] = counts.get(cluster, 0) + 1
        element_id = f"{cluster}_origin{counts[cluster]}"
        for name in halves:
            assigned[name] = element_id

    return assigned


def origin_columns(
    table: pd.DataFrame,
    contig_lengths: Dict[str, int],
) -> pd.DataFrame:
    """The ``wraps_origin`` and ``element_id`` columns for ``table``.

    ``wraps_origin`` is ``yes``/``no`` rather than ``yes``/empty: every row of a
    table that went through this module has an answer, and empty is reserved for
    the tables that never did -- a legacy four-column file, or a hand-written one
    that left the column out.
    """
    spans = spans_from_table(table)
    assigned = origin_spanning_elements(spans, contig_lengths)
    element_ids = [assigned.get(span.name, "") for span in spans]
    return pd.DataFrame(
        {
            "wraps_origin": ["yes" if element_id else "no" for element_id in element_ids],
            "element_id": element_ids,
        },
        index=table.index,
    )


def annotate(
    table: pd.DataFrame,
    reference: Path,
) -> pd.DataFrame:
    """The two columns, with the contig lengths read from ``reference``.

    This is the whole module from a back-end's point of view: hand it a clustered
    table and the genome it was called against, get back the columns.
    """
    return origin_columns(table, contig_lengths(reference))


def contig_lengths(reference: Path) -> Dict[str, int]:
    """``{contig: length}`` for every sequence in a reference FASTA."""
    with pysam.FastaFile(str(reference)) as fasta:
        return {contig: fasta.get_reference_length(contig) for contig in fasta.references}


def spans_from_table(table: pd.DataFrame) -> List[Span]:
    """The boundary-detection view of an IS table."""
    return [
        Span(
            name=str(row.is_name),
            contig=str(row.contig),
            start=int(row.start),
            stop=int(row.stop),
            cluster=str(row.cluster),
        )
        for row in table.itertuples(index=False)
    ]


def _by_contig_and_cluster(spans: Sequence[Span]) -> Dict[Tuple[str, str], List[Span]]:
    """Spans grouped by the contig and cluster they share.

    Rows with no cluster are dropped. An empty cluster column says nothing about
    which rows are one element; reading it as saying they all are would flag every
    unclustered locus near a boundary against a stranger at the far end.
    """
    grouped: Dict[Tuple[str, str], List[Span]] = {}
    for span in spans:
        if not span.cluster:
            continue
        grouped.setdefault((span.contig, span.cluster), []).append(span)
    return grouped
