"""Fake pysam-like objects for tests/test_clipped_read_search.py.

Broader than tests/fake_alignment.py's FakeAlignment (which only stands in
for ISClipped.__init__'s use of aln.references/lengths): the clipped-read
search cluster reads CIGAR strings, reference positions, and aligned pairs
off individual reads, and fetches reads per contig. These fakes cover just
that surface -- no real BAM file needed to characterize
clipped_read_search.search.
"""

from typing import List, Optional

# What ``pysam``'s ``get_reference_positions(full_length=True)`` returns: the
# reference position of every aligned base of a read, and ``None`` in place of
# every clipped one. Named so the fakes below can say it once -- a list literal
# that starts with the clipped end otherwise infers as a list of ``None`` and
# rejects the positions appended to it.
ReferencePositions = List[Optional[int]]


class FakeRead:
    def __init__(
        self,
        query_name,
        cigarstring,
        ref_positions: ReferencePositions,
        query_sequence,
        is_unmapped=False,
        is_reverse=False,
        reference_name="contig_1",
        aligned_pairs=None,
        infer_len=None,
    ):
        self.query_name = query_name
        self.cigarstring = cigarstring
        self._ref_positions = ref_positions
        self.query_sequence = query_sequence
        self.is_unmapped = is_unmapped
        self.is_reverse = is_reverse
        self.reference_name = reference_name
        self.aligned_pairs = (
            aligned_pairs
            if aligned_pairs is not None
            else [(i, p) for i, p in enumerate(ref_positions)]
        )
        self._infer_len = infer_len if infer_len is not None else len(query_sequence)

    def get_reference_positions(self, full_length=False):
        return list(self._ref_positions)

    def infer_read_length(self):
        return self._infer_len


class FakeAlignmentFetch:
    """references/lengths like FakeAlignment, plus a fetch() that returns
    canned reads per chrom regardless of start/stop -- the boundaries under
    test are chosen by the caller, not derived from real coordinate overlap.
    """

    references = ("contig_1", "contig_2")
    lengths = (2000, 2000)

    def __init__(self, reads_by_chrom):
        self._reads_by_chrom = reads_by_chrom

    def fetch(self, chrom, start, stop):
        return self._reads_by_chrom.get(chrom, [])
