#!/usr/bin/env python3
"""Rebuild ``isfinder_db.fna``, the stand-in ISFinder database.

Run from the repo root: ``python tests/fixtures/isfinder/make_isfinder_db_fixture.py``

The real ISFinder database is a few thousand sequences and is not ours to
redistribute, so the migration tests search a stand-in built from the loci the
golden table already records. Each distinct base IS name becomes one entry, named
the way real ISFinder subject ids are (``name_family_group``) with the family and
group the golden carries, and carrying the sequence of that name's longest locus.

**What it can stand in for**: extracting a locus, searching it against a
database, picking the best hit by bitscore, and splitting a subject id into name,
family and group. Those are the steps the migration back-end owns.

**What it cannot**: percent identity. A locus searched against a database whose
entry *is* that locus matches itself at 100%, where the real database -- holding
the type element rather than this genome's copy of it -- gives the golden's
98.556, 94.278 and so on. The migration tests therefore assert family, group,
cluster and the origin-spanning flags against the golden, and not ``pident``.
Reproducing ``pident`` needs the real database.
"""

import gzip
import pathlib
import sys

import pysam

sys.path.insert(0, "src")

from ijump import is_clustering, is_table  # noqa: E402

FIXTURES = pathlib.Path(__file__).parent
GOLDEN = FIXTURES.parent.parent / "goldens" / "isfinder_db_parse" / "ISTable_processing.txt"
LINE_WIDTH = 60


def main():
    golden = is_table.read_is_table(GOLDEN)

    # The longest locus of each element: a fragment would make a database entry
    # its own siblings could not match.
    longest = {}
    for row in golden.itertuples(index=False):
        base = is_clustering.base_is_name(row.is_name)
        length = int(row.stop) - int(row.start) + 1
        if base not in longest or length > longest[base][0]:
            longest[base] = (length, row)

    reference = FIXTURES / "reference.fna"
    with gzip.open(FIXTURES / "reference.fna.gz", "rb") as packed:
        reference.write_bytes(packed.read())
    pysam.faidx(str(reference))

    try:
        with pysam.FastaFile(str(reference)) as fasta:
            with open(FIXTURES / "isfinder_db.fna", "w") as out:
                for base, (_, row) in longest.items():
                    sequence = fasta.fetch(row.contig, int(row.start) - 1, int(row.stop))
                    wrapped = "\n".join(
                        sequence[i : i + LINE_WIDTH] for i in range(0, len(sequence), LINE_WIDTH)
                    )
                    out.write(f">{base}_{row.family}_{row.group}\n{wrapped}\n")
    finally:
        reference.unlink()
        (FIXTURES / "reference.fna.fai").unlink()

    print(f"wrote {len(longest)} entries to {FIXTURES / 'isfinder_db.fna'}")


if __name__ == "__main__":
    main()
