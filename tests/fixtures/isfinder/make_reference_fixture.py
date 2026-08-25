#!/usr/bin/env python3
"""Rebuild ``reference.fna.gz`` from the real assembly (isfinder-annotation 04).

    python tests/fixtures/isfinder/make_reference_fixture.py [path/to/assembly.fna]

Clustering reads the called loci out of the reference, so the parser golden needs
a reference that carries them -- but ``Test/A_baumannii_assembly.fna`` is 4 MB and
gitignored. This writes a stand-in: the same contigs at their real names and
lengths, the called loci at their real coordinates and sequences, and every other
base masked to ``N``. Locus extraction and the all-vs-all search see exactly what
they would see on the real assembly; the mask compresses to ~20 KB.

Coordinates come from the parser golden itself, so re-pinning the golden and
re-running this keeps the two in step.
"""

import os
import sys

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
DEFAULT_ASSEMBLY = os.path.join(REPO_ROOT, "Test", "A_baumannii_assembly.fna")
GOLDEN_TABLE = os.path.join(
    REPO_ROOT, "tests", "goldens", "isfinder_db_parse", "ISTable_processing.txt"
)
LINE_WIDTH = 70


def main(argv):
    assembly = argv[1] if len(argv) > 1 else DEFAULT_ASSEMBLY
    if not os.path.isfile(assembly):
        raise SystemExit(f"assembly not found: {assembly}")

    table = pd.read_csv(GOLDEN_TABLE, sep="\t")
    loci = {}
    for row in table.itertuples():
        loci.setdefault(str(row.contig), []).append((int(row.start), int(row.stop)))

    plain = os.path.join(HERE, "reference.fna")
    with pysam.FastaFile(assembly) as fasta, open(plain, "w") as out:
        for contig in fasta.references:
            if contig not in loci:
                continue
            length = fasta.get_reference_length(contig)
            masked = bytearray(b"N" * length)
            sequence = fasta.fetch(contig)
            for start, stop in loci[contig]:
                masked[start - 1 : stop] = sequence[start - 1 : stop].encode()
            out.write(f">{contig}\n")
            for offset in range(0, length, LINE_WIDTH):
                out.write(masked[offset : offset + LINE_WIDTH].decode() + "\n")

    compressed = plain + ".gz"
    pysam.tabix_compress(plain, compressed, force=True)
    os.remove(plain)
    pysam.faidx(compressed)
    print(f"wrote {compressed} (+ .fai, .gzi)")


if __name__ == "__main__":
    main(sys.argv)
