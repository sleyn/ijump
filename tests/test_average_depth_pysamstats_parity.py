"""Ticket 16: does ``ISClipped.average_depth``'s pure-``pysam`` CIGAR/span
accumulator (``isclipped.py``) reproduce ``pysamstats.load_coverage``'s true
float mean, without ``pysamstats`` installed?

Every test in this module that needs ``pysamstats`` for the reference value
calls ``pytest.importorskip("pysamstats")`` itself (not at module level), so
this file collects and passes in a plain venv with ``pysamstats`` not
installed (those tests are skipped, not failed) -- this is the actual
verification the ticket asks for -- while still giving full parity coverage
in an environment that does have ``pysamstats`` (e.g. a conda env built from
``environment.yml``).

## Round 1 (2026-08-16, closed) -- kept for the record

Tried swapping in ``pysam``'s ``count_coverage``. Filter semantics were
reconcilable (``quality_threshold=0`` closes the one real gap: pysamstats
never filters by base quality, ``count_coverage`` does by default) -- see
``test_round1_tuned_count_coverage_matches_pysamstats_filters`` and
``test_round1_tuned_count_coverage_matches_pysamstats_across_realistic_windows``.
But the round-1 bar was byte-for-byte equality with *current production
output*, and production had a pre-existing bug: ``average_depth`` called
``statistics.mean()`` directly on pysamstats' ``numpy.int32`` ``reads_all``
array, which coerces the result back to int32, truncating any fractional
mean (``mean(np.array([1, 2, 4], dtype='int32'))`` returns ``2``, not
``2.333...``). No pysam-only replacement reproduces that truncation without
deliberately reimplementing the bug, so round 1 stopped there.
``count_coverage`` also measured ~1.85x-2.11x slower than pysamstats
(recorded in the ticket's Comments, not re-measured here).

## Round 2 (2026-08-17) -- what actually shipped

A pure ``pysam`` implementation with no pileup at all: per read, clip
``[reference_start, reference_end)`` to the query window and sum the
overlap; a read's reference span is contiguous coverage for pysamstats'
purposes because every ref-consuming CIGAR op (M/=/X but also D/N) sits
inside it, and only S/I/H (which don't consume the reference) fall outside
it -- confirmed directly (``test_internal_deletion_counts_as_covered``,
``test_ref_skip_counts_as_covered``) rather than assumed. Two more
semantics had to be verified against real pysamstats output before the
swap, both confirmed by direct measurement in an env with pysamstats
installed:

- **Denominator.** ``pysamstats.load_coverage(..., pad=False)`` (the
  default, and what production used) emits a row only for reference
  positions some qualifying read actually covers -- not one for every
  position in ``[start, stop)``. The accumulator's denominator is "count of
  covered positions in the window", not window length -- see
  ``test_denominator_is_covered_positions_not_window_length``.
- **Supplementary reads.** htslib's default pileup flag filter (which
  pysamstats' pileup is built on) excludes UNMAP/SECONDARY/QCFAIL/DUP but
  *not* SUPPLEMENTARY -- confirmed empirically, not assumed from the flag
  mask alone. iJump is a transposon-insertion caller and supplementary
  alignments of clipped reads are precisely its core signal, so
  ``average_depth`` keeps them -- see
  ``test_supplementary_reads_are_counted``.

This closes both round-1 blockers at the root: no base-quality filter to
reconcile (the arithmetic never reads ``SEQ``/qualities), and no per-column
pileup to be slow (a CIGAR-span sum measured *faster* than pysamstats, not
just under the ~1.5x bar -- see the ticket's Comments for numbers).

Fixing the truncation was in scope for round 2 (it's what any correct
replacement does incidentally) -- see
``test_production_average_depth_no_longer_truncates_fractional_means``,
which inverts round 1's finding now that it's fixed.
"""

import random
from statistics import StatisticsError, mean

import pysam
import pytest

from ijump.isclipped import ISClipped

HEADER_TEMPLATE = {"HD": {"VN": "1.6", "SO": "coordinate"}}


def _make_read(name, pos, cigar, flag, qual=40, ref_id=0):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigar = cigar
    qlen = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    a.query_sequence = "A" * qlen
    a.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * qlen)
    return a


def _simple_read(name, pos, seqlen, qual, flag, ref_id=0):
    return _make_read(name, pos, [(0, seqlen)], flag, qual=qual, ref_id=ref_id)


def _build_bam(path, contig, contig_len, reads):
    header = dict(HEADER_TEMPLATE, SQ=[{"SN": contig, "LN": contig_len}])
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in sorted(reads, key=lambda r: r.reference_start):
            out.write(r)
    pysam.index(str(path))


def _count_coverage_average(aln, chrom, start, stop, **kwargs):
    aln_depth = aln.count_coverage(chrom, start, stop, **kwargs)
    depth = sum(map(sum, aln_depth))
    return depth / len(aln_depth[0])


def _pysamstats_true_mean(pysamstats, aln, chrom, start, stop):
    """True (unrounded) mean of pysamstats' reads_all, or None if the
    window has zero covered positions (pysamstats/pad=False emits no rows
    at all in that case)."""
    c = pysamstats.load_coverage(
        aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000
    )
    if len(c.reads_all) == 0:
        return None
    return sum(int(x) for x in c.reads_all) / len(c.reads_all)


# --- Round 1 (historical, kept for the record) -----------------------------


def test_round1_tuned_count_coverage_matches_pysamstats_filters(tmp_path):
    """Region has: two normal reads, one low-base-quality read (qual=5,
    below count_coverage's default quality_threshold=15), and one
    duplicate-flagged read. True per-position depth (over reads that
    survive filtering) is uniform, so the mean is an integer -- isolating
    filter behavior from the truncation issue covered separately below.
    """
    pysamstats = pytest.importorskip("pysamstats")

    bam = tmp_path / "filters.bam"
    contig, contig_len = "chr1", 200
    reads = [
        _simple_read("r1_normal", 10, 50, 40, 0),
        _simple_read("r2_normal", 10, 50, 40, 0),
        _simple_read("r3_lowbaseq", 10, 50, 5, 0),
        _simple_read("r4_duplicate", 10, 50, 40, 1024),
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    chrom, start, stop = contig, 20, 40

    c = pysamstats.load_coverage(
        aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000
    )
    pysamstats_avg = mean(c.reads_all)
    # r4 is excluded (duplicate-flagged); r3 (low base quality) IS counted
    # -- pysamstats never filters by base quality. r1, r2, r3 -> depth 3.
    assert pysamstats_avg == 3

    # count_coverage's own default quality_threshold=15 disagrees: it also
    # drops r3 for low base quality, leaving only r1, r2 -> depth 2.
    default_avg = _count_coverage_average(aln, chrom, start, stop)
    assert default_avg == 2
    assert default_avg != pysamstats_avg

    # Tuned: quality_threshold=0 restores r3.
    tuned_avg = _count_coverage_average(aln, chrom, start, stop, quality_threshold=0)
    assert tuned_avg == pysamstats_avg == 3


def test_round1_tuned_count_coverage_matches_pysamstats_across_realistic_windows(tmp_path):
    """~80x coverage over a 5kb window, many overlapping reads. Confirms
    filter-level parity holds beyond the single hand-built region above."""
    pysamstats = pytest.importorskip("pysamstats")

    random.seed(11)
    bam = tmp_path / "sweep.bam"
    contig, contig_len = "chr1", 5000
    read_len = 100
    n_reads = int(contig_len * 80 / read_len)
    reads = [
        _simple_read(f"r{i}", random.randint(0, contig_len - read_len), read_len, 40, 0)
        for i in range(n_reads)
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))

    mismatches = 0
    for _ in range(50):
        w = random.choice((50, 100, 250))
        start = random.randint(0, contig_len - w - 1)
        stop = start + w

        pysamstats_true_mean = _pysamstats_true_mean(pysamstats, aln, contig, start, stop)
        tuned_avg = _count_coverage_average(aln, contig, start, stop, quality_threshold=0)
        if pysamstats_true_mean != tuned_avg:
            mismatches += 1

    assert mismatches == 0


def test_production_average_depth_no_longer_truncates_fractional_means(tmp_path):
    """Round 1's stop condition, inverted: the truncation bug that made
    round 1's numerically-correct candidates disagree with production is
    now fixed as part of round 2's swap (ticket 16 Step 4). Production
    returns the true 1.5, not the old int32-coerced 1.
    """
    bam = tmp_path / "fractional.bam"
    contig, contig_len = "chr1", 200
    # Read A spans [10, 20): depth 1 across all 10 positions.
    # Read B spans [10, 15): depth +1 across the first 5 positions only.
    reads = [
        _simple_read("A", 10, 10, 40, 0),
        _simple_read("B", 10, 5, 40, 0),
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    chrom, start, stop = contig, 10, 20

    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")
    production_avg = isc.average_depth(chrom, start, stop)

    assert production_avg == 1.5
    assert isinstance(production_avg, float)


# --- Round 2: new semantics against real production average_depth ---------


def test_supplementary_reads_are_counted(tmp_path):
    """htslib's default pileup flag filter (pysamstats' basis) excludes
    UNMAP/SECONDARY/QCFAIL/DUP but not SUPPLEMENTARY. Clipped supplementary
    alignments are iJump's core signal, so they must stay counted.
    """
    bam = tmp_path / "supplementary.bam"
    contig, contig_len = "chr1", 200
    reads = [
        _simple_read("r1", 10, 20, 40, 0),
        _simple_read("r2_supplementary", 10, 20, 40, 2048),
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    assert isc.average_depth(contig, 15, 25) == 2.0


def test_internal_deletion_counts_as_covered(tmp_path):
    """A read spanning a CIGAR 'D' still occupies those reference
    positions in pysamstats' reads_all (is_del pileup rows still count the
    read). A read.reference_start:reference_end span naturally includes
    the deleted region, so this must match without special-casing D.
    """
    bam = tmp_path / "deletion.bam"
    contig, contig_len = "chr1", 200
    # 5M 5D 5M: ref span [10, 25), deletion at [15, 20).
    reads = [_make_read("d1", 10, [(0, 5), (2, 5), (0, 5)], 0)]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    assert isc.average_depth(contig, 10, 25) == 1.0


def test_ref_skip_counts_as_covered(tmp_path):
    """Same as the deletion case, for CIGAR 'N' (reference skip)."""
    bam = tmp_path / "refskip.bam"
    contig, contig_len = "chr1", 200
    reads = [_make_read("n1", 10, [(0, 5), (3, 5), (0, 5)], 0)]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    assert isc.average_depth(contig, 10, 25) == 1.0


def test_denominator_is_covered_positions_not_window_length(tmp_path):
    """pysamstats/pad=False emits a row only for positions that are
    actually covered by some read, not one per position in the window --
    average_depth's mean must be over covered positions, not window
    length, or a window with partial zero-coverage flanks would silently
    dilute the mean.
    """
    bam = tmp_path / "partial.bam"
    contig, contig_len = "chr1", 200
    # Read covers only [10, 20); window queried is [0, 30), which has 20
    # zero-coverage positions flanking it.
    reads = [_simple_read("r1", 10, 10, 40, 0)]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    # Mean over the 10 covered positions (depth 1 throughout) is 1.0, not
    # 10/30 == 0.333... (mean over the full window length).
    assert isc.average_depth(contig, 0, 30) == 1.0


def test_zero_coverage_window_raises(tmp_path):
    """A window with no covered positions at all has no true mean --
    pysamstats' own reads_all array would be empty there too (mean() of an
    empty pysamstats array already raised StatisticsError in production
    before this ticket; the pysam-only accumulator preserves that)."""
    bam = tmp_path / "empty.bam"
    contig, contig_len = "chr1", 200
    reads = [_simple_read("r1", 10, 10, 40, 0)]  # covers [10, 20) only
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    with pytest.raises(StatisticsError):
        isc.average_depth(contig, 100, 110)


def _build_mixed_feature_bam(path, seed):
    """~80x coverage over a 200kb contig, with internal deletions,
    ref-skips, and supplementary/duplicate/secondary flags mixed into
    ~12% of reads -- exercises every round-2 semantics axis at once."""
    random.seed(seed)
    contig, contig_len = "chr1", 200000
    read_len = 100
    n_reads = int(contig_len * 80 / read_len)
    reads = []
    for i in range(n_reads):
        pos = random.randint(0, contig_len - read_len - 20)
        flag = 0
        r = random.random()
        if r < 0.03:
            flag = 2048  # supplementary
        elif r < 0.05:
            flag = 1024  # duplicate
        elif r < 0.06:
            flag = 256  # secondary
        if flag == 0 and r < 0.10:
            gap = random.randint(1, 20)
            cigar = [(0, 40), (2, gap), (0, read_len - 40)]  # internal deletion
        elif flag == 0 and r < 0.12:
            gap = random.randint(1, 20)
            cigar = [(0, 40), (3, gap), (0, read_len - 40)]  # ref skip
        else:
            cigar = [(0, read_len)]
        reads.append(_make_read(f"r{i}", pos, cigar, flag))
    _build_bam(path, contig, contig_len, reads)
    return contig, contig_len


@pytest.mark.parametrize("n_windows", [50, 300])
def test_production_matches_pysamstats_across_realistic_windows_with_gaps_and_supplementary(
    tmp_path, n_windows
):
    """Sweep against the real production ISClipped.average_depth, at both
    scales round 1 used (50-region and 300-region sweeps; ticket 16 Step
    2's bar). 0 mismatches required against pysamstats' true float mean,
    across windows that include zero-coverage stretches."""
    pysamstats = pytest.importorskip("pysamstats")

    bam = tmp_path / "sweep.bam"
    contig, contig_len = _build_mixed_feature_bam(bam, seed=42)
    aln = pysam.AlignmentFile(str(bam))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    mismatches = 0
    for _ in range(n_windows):
        w = random.choice((50, 100, 250, 500))
        start = random.randint(0, contig_len - w - 1)
        stop = start + w

        pysamstats_true_mean = _pysamstats_true_mean(pysamstats, aln, contig, start, stop)
        try:
            production_avg = isc.average_depth(contig, start, stop)
        except StatisticsError:
            production_avg = None

        if pysamstats_true_mean != production_avg:
            mismatches += 1

    assert mismatches == 0
