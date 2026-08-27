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

import hashlib
import os
from typing import Dict, Tuple, Union

import pandas as pd

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# Column order of the table as written. The first four are the legacy columns and
# keep their positions, so a legacy table is the new table minus its tail.
COLUMNS = (
    "is_name",
    "contig",
    "start",
    "stop",
    "family",
    "group",
    "cluster",
    "pident",
    "wraps_origin",
    "element_id",
)

# Columns a legacy headerless table carries, in order.
LEGACY_COLUMNS = COLUMNS[:4]

# How an operator turns a table with no cluster assignment into one that has it.
# Named in the error below rather than left for the reader to find, because that
# error is the only thing standing between a legacy table and a working run.
MIGRATE_SUBCOMMAND = "ijump migrate-is-table"


# How a report names the IS table it was built from. Short enough to read in a
# file header, long enough that two different tables will not collide.
FINGERPRINT_LENGTH = 16


class MissingClusterColumn(Exception):
    """An IS table carries no cluster where the caller requires one.

    Raised rather than filled in: the cluster is the answer to which loci a
    clipped read cannot tell apart, and nothing in the table's other columns
    answers that. Guessing it from the IS name is exactly what
    isfinder-annotation 06 removed.
    """


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


def with_all_columns(table: pd.DataFrame) -> pd.DataFrame:
    """``table`` widened to all of COLUMNS, in order, absent ones empty.

    The one place the fill rule lives: a back-end states the columns it can
    speak to and this fills in the rest, rather than each back-end carrying its
    own list of empty strings to fall out of step with COLUMNS. Columns beyond
    COLUMNS are kept, at the end -- a table may carry more than its reader knows
    about, and dropping such a column on the way through would defeat the header.
    """
    widened = table.copy()
    for column in COLUMNS:
        if column not in widened.columns:
            widened[column] = ""
    extra = [column for column in widened.columns if column not in COLUMNS]
    return widened[list(COLUMNS) + extra]


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
        return with_all_columns(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))

    # Legacy tables were split on arbitrary whitespace, so keep accepting that.
    table = pd.read_csv(
        path,
        sep=r"\s+",
        names=list(LEGACY_COLUMNS),
        dtype=str,
        keep_default_na=False,
    )
    return with_all_columns(table)


def cluster_by_name(table: pd.DataFrame) -> Dict[str, str]:
    """Map each row's ``is_name`` to its ``cluster``.

    The lookup callers group on: two rows sharing a cluster are copies or
    fragments of one mobile element, whatever their names say. A row with the
    cluster left empty -- every row of a legacy table, or one cell an operator
    blanked by hand -- raises rather than being grouped on its name.

    Two rows sharing an ``is_name`` raise too, for the reason
    ``is_clustering.loci_from_table`` rejects them: a dict keyed on the name
    would keep the last row's cluster and drop the other silently.
    """
    names = [str(row.is_name) for row in table.itertuples(index=False)]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(
            "IS table has more than one row named: " + ", ".join(duplicates) + ". "
            "Element names have to be unique."
        )

    missing = [
        str(row.is_name) for row in table.itertuples(index=False) if not str(row.cluster).strip()
    ]
    if missing:
        raise MissingClusterColumn(
            "The IS table has no cluster assigned for: "
            + ", ".join(missing)
            + ". The cluster column says which loci are copies of one element, and "
            "grouping cannot be inferred from the IS names. Produce a table that "
            f"carries it with `{MIGRATE_SUBCOMMAND}`, or fill the column in by hand."
        )

    return {str(row.is_name): str(row.cluster) for row in table.itertuples(index=False)}


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


def fingerprint(table: pd.DataFrame) -> str:
    """A short digest identifying which IS table a run was annotated against.

    Cluster names are *derived* from the loci, not fixed labels, so the same name
    can mean different elements in two runs annotated against different tables or
    different references. Reports carry this digest so a multi-sample merge can
    tell that its samples share one vocabulary before it joins them on it
    (isfinder-annotation 07).

    Taken over the canonical columns only, so that reading a legacy table into
    the wider shape -- which fills the annotation columns with empty strings --
    gives the same answer as reading a headered table that omits them. Extra
    columns a reader does not know about are left out for the same reason: they
    do not change what the cluster names mean.
    """
    canonical = table[list(COLUMNS)].to_csv(sep="\t", header=True, index=False)
    return hashlib.sha256(canonical.encode()).hexdigest()[:FINGERPRINT_LENGTH]
