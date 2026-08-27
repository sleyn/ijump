"""A minimal stand-in for pysam.AlignmentFile.

ISClipped.__init__ (isclipped.py:21-119) only touches ``aln.references`` and
``aln.lengths`` (isclipped.py:80-83), so constructing an ISClipped instance in
a test needs no real BAM file.
"""


class FakeAlignment:
    references = ("contig_1",)
    lengths = (10000,)
