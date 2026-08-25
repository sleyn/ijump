"""Read and write the IS table -- the IS elements a run works from, one row each.

``ijump isfinder-db-parse`` writes it, ``ijump run`` reads it, and an operator may
edit it in between; it is the one file both halves of the codebase agree on, so
its format lives here rather than being spelled out at each end.

The format is a headered TSV. A header names each column, so a table can be read
without knowing which generation wrote it and later columns can be added without
breaking readers -- which matters, since the annotation columns below are the
first of several. Tables written before the header existed (four
whitespace-separated columns: is_name, contig, start, stop) are recognised by the
absence of the header row and read with the added columns empty.

Every field is kept as a string. Coordinates are consumed as text in some places
and as ``int`` in others, and text is the form that survives a round trip through
the file unchanged.
"""

import os
from typing import Tuple, Union

import pandas as pd

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# Column order of the table as written. The first four are the legacy columns and
# keep their positions, so a legacy table is the new table minus its tail.
COLUMNS = ("is_name", "contig", "start", "stop", "family", "group", "pident")

# Columns a legacy headerless table carries, in order.
LEGACY_COLUMNS = COLUMNS[:4]


def parse_subject_id(sseqid: str) -> Tuple[str, str, str]:
    """Split an ISFinder BLAST subject id into ``(name, family, group)``.

    Subject ids look like ``ISAba18_IS3_IS51`` -- family and group are the last
    two underscore-separated fields, so the split runs from the right. Eleven
    database entries have an underscore inside the *name* itself
    (``ISBj2_B_IS5_IS5``); splitting from the left truncated those.

    An id with fewer than three fields carries no annotation to recover, so it is
    kept whole as the name.
    """
    parts = sseqid.rsplit("_", 2)
    if len(parts) < 3:
        return sseqid, "", ""
    name, family, group = parts
    return name, family, group


def write_is_table(table: pd.DataFrame, path: Path) -> None:
    """Write ``table`` as the headered TSV, in the canonical column order."""
    table.to_csv(path, sep="\t", columns=list(COLUMNS), header=True, index=False)


def read_is_table(path: Path) -> pd.DataFrame:
    """Read an IS table, headered or legacy, into a frame carrying all of COLUMNS.

    Columns the file does not carry come back as empty strings -- a legacy table
    has none of the annotation ones, and a hand-written table may leave out any
    it has nothing to say about. Columns beyond COLUMNS are kept as they were
    read: a header exists so that a table can carry more than its reader knows
    about, and silently dropping such a column on the way through would defeat
    that.
    """
    if _has_header(path):
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for column in COLUMNS:
            if column not in table.columns:
                table[column] = ""
        extra = [column for column in table.columns if column not in COLUMNS]
        return table[list(COLUMNS) + extra]

    # Legacy tables were split on arbitrary whitespace, so keep accepting that.
    table = pd.read_csv(
        path,
        sep=r"\s+",
        names=list(LEGACY_COLUMNS),
        dtype=str,
        keep_default_na=False,
    )
    for column in COLUMNS[len(LEGACY_COLUMNS) :]:
        table[column] = ""
    return table[list(COLUMNS)]


def _has_header(path: Path) -> bool:
    """Whether the file starts with the header row.

    Matched on the whole leading run of column names rather than on the first
    cell alone, so a legacy table whose first row happens to be named
    ``is_name`` is still read as data.
    """
    with open(path, "r") as table_file:
        first_line = table_file.readline()
    fields = first_line.rstrip("\n").split("\t")
    return fields[: len(LEGACY_COLUMNS)] == list(LEGACY_COLUMNS)
