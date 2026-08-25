"""Group the called IS loci into clusters -- the copies of one mobile element.

The IS table's ``is_name`` is the nearest ISFinder database entry plus a copy
number the parser invents, and neither half is a reliable way to tell which rows
are the same element. Two fragments of one element land on different database
entries and get different names (``IS17_1``, ``IS17_2`` and ``ISAba12_1`` are one
ISAba12-like copy with the middle 819 bp absent from the assembly), while two
elements 83% apart share a name whenever the database has no closer entry for
one of them. So the grouping is computed from the sequences themselves: extract
each locus from the reference, align them all against each other, and link the
pairs that a clipped read could not tell apart.

**Coverage is measured on the shorter locus.** A 76 bp remnant aligning over its
whole length to a 1039 bp element is the same element as far as this pipeline can
tell -- a read clipped at the remnant is indistinguishable from one clipped at the
full copy -- so the remnant joins, even though the alignment covers 7% of its
parent.

**Linkage is single, not complete.** ``IS17_1`` and ``IS17_2`` are opposite ends
of one element and do not align to each other at all; they reach each other only
through ``ISAba12_1``. Complete linkage would break the exact case this exists to
fix. The price is chaining -- a fragment landing in a stretch conserved between two
distinct elements would merge them -- so every internal pair that fails the
threshold is logged by name. Nothing in the alignment distinguishes a wanted chain
from an unwanted one; the operator reads the log and edits the cluster column
before running the pipeline.
"""

import argparse
import logging
import os
import re
import string
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Set, Tuple, Union

import pandas as pd
import pysam

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# The blastn flags every locus search shares.
#
# The sensitive ``blastn`` task rather than the binary's megablast default: the
# shortest loci here are under 100 bp, and megablast's seed of 28 would need an
# exact 28-mer to find them at all. Extra weak hits are harmless -- callers apply
# their own thresholds. One list rather than one per search, because this is one
# tuning decision about one kind of query, and two copies would drift.
SENSITIVE_BLASTN = (
    "-task",
    "blastn",
    "-word_size",
    "11",
    "-dust",
    "no",
    "-evalue",
    "1e-5",
)

# Minimum percent identity of an alignment for it to link two loci.
IDENTITY_DEFAULT = 95.0
# Minimum fraction of the *shorter* locus the alignment has to span.
COVERAGE_DEFAULT = 0.8

# The copy number ``isfinder_db_parcer`` appends to a database name: strictly a
# trailing underscore and digits, so the 11 database names carrying an underscore
# of their own (``ISBj2_B``) keep it.
_COPY_SUFFIX = re.compile(r"_\d+$")


def threshold_type(low: float, high: float, unit: str):
    """An argparse type for a cluster threshold, rejecting values off its scale.

    The two thresholds are on different scales -- identity is a percent, coverage
    a fraction, because that is how BLAST reports each -- and an out-of-scale
    value does not fail, it just silently answers wrongly: ``--cluster-coverage
    80`` would leave every element in a cluster of its own. Shared by every
    back-end's CLI so they cannot disagree on what a threshold is.
    """

    def parse(value):
        number = float(value)
        if not low <= number <= high:
            raise argparse.ArgumentTypeError(
                f"{value} is not {unit}; expected a value between {low} and {high}"
            )
        return number

    return parse


@dataclass(frozen=True)
class Locus:
    """One row of the IS table, as clustering sees it."""

    name: str
    contig: str
    start: int
    length: int


@dataclass(frozen=True)
class Hit:
    """One alignment between two loci: a single BLAST HSP.

    HSPs are not summed. A cluster forms on evidence that a clipped read of any
    reasonable length could match either locus, and that is one contiguous
    alignment, not several scattered ones.
    """

    query: str
    subject: str
    identity: float
    aligned_length: int


@dataclass
class Cluster:
    """Loci that are the same mobile element.

    ``members`` is ordered longest first, so ``members[0]`` is the representative
    the cluster is named after. ``chained_pairs`` are the internal pairs that do
    not meet the threshold themselves -- empty for a cluster every pair of which
    aligns.
    """

    members: List[Locus]
    chained_pairs: List[Tuple[str, str]] = field(default_factory=list)

    @property
    def representative(self) -> Locus:
        return self.members[0]


def cluster_column(
    loci: Sequence[Locus],
    hits: Iterable[Hit],
    identity: float = IDENTITY_DEFAULT,
    coverage: float = COVERAGE_DEFAULT,
) -> Dict[str, str]:
    """Cluster ``loci`` and return ``{locus name: cluster name}``.

    Chained clusters are logged here rather than at the call site, so every
    back-end reaching this function reports them the same way.
    """
    clusters = cluster_loci(loci, hits, identity, coverage)
    names = name_clusters(clusters)

    for cluster, name in zip(clusters, names):
        for first, second in cluster.chained_pairs:
            logging.warning(
                "Cluster %s is chained: %s and %s share it without meeting "
                "%.4g%% identity over %.4g%% of the shorter locus. Edit the "
                "cluster column if they are different elements.",
                name,
                first,
                second,
                identity,
                coverage * 100,
            )

    return {locus.name: name for cluster, name in zip(clusters, names) for locus in cluster.members}


def cluster_loci(
    loci: Sequence[Locus],
    hits: Iterable[Hit],
    identity: float = IDENTITY_DEFAULT,
    coverage: float = COVERAGE_DEFAULT,
) -> List[Cluster]:
    """Single-linkage clusters of ``loci`` over the alignments in ``hits``.

    Pure: every alignment it works from is passed in, so the linkage rule is
    testable against the measured numbers without running BLAST.
    """
    by_name = {locus.name: locus for locus in loci}
    linked = _linked_pairs(by_name, hits, identity, coverage)

    parent = {name: name for name in by_name}
    for first, second in linked:
        _union(parent, first, second)

    grouped: Dict[str, List[Locus]] = {}
    for locus in loci:
        grouped.setdefault(_find(parent, locus.name), []).append(locus)

    clusters = [Cluster(members=_by_representative(members)) for members in grouped.values()]
    for cluster in clusters:
        cluster.chained_pairs = _chained_pairs(cluster, linked)

    # Order clusters the way they are named, so the result does not depend on
    # dictionary insertion order.
    clusters.sort(key=lambda cluster: _coordinate_key(cluster.representative))
    return clusters


def name_clusters(clusters: Sequence[Cluster]) -> List[str]:
    """Name each cluster after its representative, in the order given.

    A cluster is named for the base IS name of its longest member -- the full
    copy, not a fragment of it. Two clusters can want the same name (the database
    is sparse relative to real IS diversity, so loci ~85% apart share a nearest
    entry), and a duplicate would be a correctness bug rather than a cosmetic one:
    ``ISClipped.is_coords`` is a dict keyed on this name, and precise mode groups
    on it. So colliding names -- and only colliding names -- take a suffix,
    largest cluster first, then coordinate.
    """
    wanted = [base_is_name(cluster.representative.name) for cluster in clusters]

    collisions = {name for name in wanted if wanted.count(name) > 1}
    ranked = sorted(
        (index for index, name in enumerate(wanted) if name in collisions),
        key=lambda index: (
            -len(clusters[index].members),
            _coordinate_key(clusters[index].representative),
        ),
    )

    named = list(wanted)
    seen: Dict[str, int] = {}
    for index in ranked:
        name = wanted[index]
        named[index] = f"{name}.{_suffix(seen.get(name, 0))}"
        seen[name] = seen.get(name, 0) + 1
    return named


def base_is_name(is_name: str) -> str:
    """The IS name without the copy number the parser appended."""
    return _COPY_SUFFIX.sub("", is_name)


def loci_from_table(table: pd.DataFrame) -> List[Locus]:
    """The clustering view of an IS table.

    Two rows sharing an ``is_name`` are rejected rather than clustered. Everything
    downstream keys on that name -- linkage here, ``ISClipped.is_coords`` later --
    so a duplicate would quietly drop one of the two rows instead of failing.
    """
    loci = [
        Locus(
            name=str(row.is_name),
            contig=str(row.contig),
            start=int(row.start),
            length=int(row.stop) - int(row.start) + 1,
        )
        for row in table.itertuples(index=False)
    ]

    names = [locus.name for locus in loci]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(
            "IS table has more than one row named: " + ", ".join(duplicates) + ". "
            "Element names have to be unique."
        )
    return loci


def extract_loci(reference: Path, loci: Sequence[Locus]) -> Dict[str, str]:
    """Read each locus's sequence out of the reference FASTA.

    Table coordinates are 1-based and inclusive, as BLAST reports them.
    """
    sequences = {}
    with pysam.FastaFile(str(reference)) as fasta:
        available = set(fasta.references)
        for locus in loci:
            if locus.contig not in available:
                raise ValueError(
                    f"contig {locus.contig!r} of IS element {locus.name!r} is not in "
                    f"{reference} -- is this the reference the IS table was called against?"
                )
            sequences[locus.name] = fasta.fetch(
                locus.contig, locus.start - 1, locus.start + locus.length - 1
            )
    return sequences


def all_vs_all_search(sequences: Dict[str, str]) -> List[Hit]:
    """Align every locus against every other with ``blastn``.

    BLAST+ is already a hard dependency of a run, so this adds nothing to
    install. Sensitivity comes from SENSITIVE_BLASTN above; the identity and
    coverage thresholds discard the extra weak hits it lets through.

    ``-max_target_seqs`` is raised to the number of loci. Its default of 500
    truncates the per-query hit list, and a truncated list here is not a slower
    answer but a wrong one: single linkage would silently lose an edge on a
    genome carrying more than 500 called elements, which is the conflation this
    whole module exists to prevent.
    """
    if len(sequences) < 2:
        return []

    with tempfile.TemporaryDirectory() as workdir:
        query = os.path.join(workdir, "loci.fna")
        write_query_fasta(sequences, query)

        database = os.path.join(workdir, "loci")
        out_file = os.path.join(workdir, "loci_vs_loci.out")
        run_blast_command(["makeblastdb", "-in", query, "-dbtype", "nucl", "-out", database])
        run_blast_command(
            [
                "blastn",
                "-query",
                query,
                "-db",
                database,
                "-out",
                out_file,
                "-outfmt",
                "6 qseqid sseqid pident length",
                "-max_target_seqs",
                str(len(sequences)),
                *SENSITIVE_BLASTN,
            ]
        )
        return _read_hits(out_file)


def annotate(
    table: pd.DataFrame,
    reference: Path,
    identity: float = IDENTITY_DEFAULT,
    coverage: float = COVERAGE_DEFAULT,
) -> pd.Series:
    """The ``cluster`` column for ``table``, computed against ``reference``.

    This is the whole subsystem from a back-end's point of view: hand it a table
    with names and coordinates, get back the column.
    """
    loci = loci_from_table(table)
    hits = all_vs_all_search(extract_loci(reference, loci))
    column = cluster_column(loci, hits, identity, coverage)
    return pd.Series([column[locus.name] for locus in loci], index=table.index)


# --- linkage ---------------------------------------------------------------


def _linked_pairs(
    by_name: Dict[str, Locus],
    hits: Iterable[Hit],
    identity: float,
    coverage: float,
) -> Set[Tuple[str, str]]:
    """The unordered pairs whose alignment meets both thresholds.

    blastn reports each pair from both ends, and the two reports need not be
    identical, so a pair links if *either* direction qualifies.
    """
    linked = set()
    for hit in hits:
        if hit.query == hit.subject:
            continue
        query, subject = by_name.get(hit.query), by_name.get(hit.subject)
        if query is None or subject is None:
            continue
        if hit.identity < identity:
            continue
        if hit.aligned_length < coverage * min(query.length, subject.length):
            continue
        linked.add(_pair(hit.query, hit.subject))
    return linked


def _pair(first: str, second: str) -> Tuple[str, str]:
    return (first, second) if first <= second else (second, first)


def _find(parent: Dict[str, str], name: str) -> str:
    while parent[name] != name:
        parent[name] = parent[parent[name]]
        name = parent[name]
    return name


def _union(parent: Dict[str, str], first: str, second: str) -> None:
    first_root, second_root = _find(parent, first), _find(parent, second)
    if first_root != second_root:
        parent[second_root] = first_root


def _chained_pairs(cluster: Cluster, linked: Set[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Internal pairs of ``cluster`` that do not meet the threshold themselves."""
    names = [locus.name for locus in cluster.members]
    return sorted(
        _pair(first, second)
        for index, first in enumerate(names)
        for second in names[index + 1 :]
        if _pair(first, second) not in linked
    )


# --- naming -----------------------------------------------------------------


def _by_representative(members: List[Locus]) -> List[Locus]:
    """Cluster members, longest first.

    Length is what picks the representative: a cluster of a full copy and its own
    fragments is that element, so the full copy names it. Identical copies leave
    length settling nothing, and the coordinate breaks the tie -- a rule rather
    than the order rows happened to arrive in.
    """
    return sorted(members, key=lambda locus: (-locus.length, locus.contig, locus.start))


def _coordinate_key(locus: Locus) -> Tuple[str, int]:
    """Where a locus sits, as a sort key -- the last tie-break in every ordering
    here, so that two runs over the same table agree on cluster names."""
    return (locus.contig, locus.start)


def _suffix(index: int) -> str:
    """``a``, ``b``, ... ``z``, ``aa``, ``ab``: bijective base-26.

    Twenty-six colliding clusters is not a realistic genome, but wrapping back
    onto an earlier letter would reintroduce the duplicate this exists to prevent.
    """
    letters = string.ascii_lowercase
    suffix = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, len(letters))
        suffix = letters[remainder] + suffix
    return suffix


# --- BLAST ------------------------------------------------------------------


def write_query_fasta(sequences: Dict[str, str], path: Path) -> None:
    """Write locus sequences as the query FASTA a blastn search reads."""
    with open(path, "w") as fasta:
        for name, sequence in sequences.items():
            fasta.write(f">{name}\n{sequence}\n")


def run_blast_command(command: List[str]) -> None:
    """Run one BLAST+ command, turning its failures into logged errors.

    Public because more than one back-end searches.
    """
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        logging.error("%s not found. Is BLAST+ installed and on PATH?", command[0])
        raise
    except subprocess.CalledProcessError as error:
        logging.error(
            "%s failed (exit %d): %s",
            command[0],
            error.returncode,
            error.stderr.decode(errors="replace"),
        )
        raise


def _read_hits(out_file: Path) -> List[Hit]:
    hits = []
    with open(out_file) as blast_out:
        for line in blast_out:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            query, subject, identity, aligned_length = fields[:4]
            hits.append(Hit(query, subject, float(identity), int(aligned_length)))
    return hits
