"""Ticket 16 characterization test: can pysam's ``count_coverage`` exactly
replace ``pysamstats.load_coverage`` inside ``ISClipped.average_depth``
(``isclipped.py:643-649``)?

Two things had to be checked, both against real ``pysam.AlignmentFile``
objects (synthetic BAMs built here, since ``tests/fixtures/tiny.bam`` is too
uniform -- 7 reads, all flag 0, all base quality 40 -- to exercise either
library's read-level filtering):

1. **Filter semantics.** ``count_coverage``'s default ``read_callback='all'``
   already matches what ``pysamstats.load_coverage`` does: both exclude
   unmapped/secondary/qcfail/duplicate-flagged reads (pysamstats' own
   ``no_dup`` toggle turns out to have no effect -- verified directly by
   feeding it a duplicate-only region and varying ``no_dup`` -- because the
   underlying ``pysam`` pileup it's built on already applies its own default
   flag filter before pysamstats' knob is consulted). The one real
   discrepancy is base-quality filtering: ``count_coverage`` filters bases
   below phred 15 by default (``quality_threshold=15``); pysamstats never
   filters by base quality. Passing ``quality_threshold=0`` to
   ``count_coverage`` (default ``read_callback`` otherwise) closes this gap
   exactly -- see ``test_tuned_count_coverage_matches_pysamstats_filters``
   and ``test_tuned_count_coverage_matches_pysamstats_across_realistic_windows``.

2. **Numeric result.** Even with filters reconciled, the two candidate
   implementations differ from what ``ISClipped.average_depth`` actually
   returns *today*, because of a pre-existing quirk: pysamstats' coverage
   arrays are ``numpy.int32``, and Python's ``statistics.mean()`` -- which
   ``average_depth`` calls directly on that array -- detects the element
   type and coerces its final result back to it. For a numpy int32 array
   this *truncates* any fractional mean instead of returning a float:
   ``statistics.mean(np.array([1, 2, 4], dtype='int32'))`` returns ``2``,
   not ``2.333...``. A straightforward ``count_coverage``-based replacement
   (plain ``depth / len(...)`` division, matching the commented-out
   pysam-only code already sitting above the current implementation) does
   not reproduce that truncation -- it returns the mathematically correct
   float. See ``test_production_average_depth_truncates_fractional_means``.

Net: filter-level parity is achievable, but the *current production output*
is not reproducible without deliberately reimplementing an int32-truncation
bug, which is out of this ticket's scope. See
``.scratch/isclipped-refactor/issues/16-evaluate-dropping-pysamstats.md``
for the full characterization and the decision not to swap.
"""

from statistics import mean

import pysam
import pysamstats

from ijump.isclipped import ISClipped

HEADER_TEMPLATE = {"HD": {"VN": "1.6", "SO": "coordinate"}}


def _make_read(name, pos, seqlen, qual, flag, ref_id=0):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "A" * seqlen
    a.flag = flag
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigar = [(0, seqlen)]
    a.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * seqlen)
    return a


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


def _pysamstats_average(aln, chrom, start, stop):
    c = pysamstats.load_coverage(
        aln, chrom=chrom, start=start, end=stop, truncate=True, max_depth=300000
    )
    return mean(c.reads_all)


# --- 1. Filter semantics: reconcilable ------------------------------------


def test_tuned_count_coverage_matches_pysamstats_filters(tmp_path):
    """Region has: two normal reads, one low-base-quality read (qual=5,
    below count_coverage's default quality_threshold=15), and one
    duplicate-flagged read. True per-position depth (over reads that
    survive filtering) is uniform, so the mean is an integer -- isolating
    filter behavior from the separate int32-truncation issue covered below.
    """
    bam = tmp_path / "filters.bam"
    contig, contig_len = "chr1", 200
    reads = [
        _make_read("r1_normal", 10, 50, 40, 0),
        _make_read("r2_normal", 10, 50, 40, 0),
        _make_read("r3_lowbaseq", 10, 50, 5, 0),
        _make_read("r4_duplicate", 10, 50, 40, 1024),
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    chrom, start, stop = contig, 20, 40

    pysamstats_avg = _pysamstats_average(aln, chrom, start, stop)
    # r4 is excluded (pysamstats' pileup already drops duplicate-flagged
    # reads regardless of its no_dup toggle -- verified separately by
    # feeding it a duplicate-only region and varying no_dup, no effect).
    # r3 (low base quality) IS counted: pysamstats never filters by base
    # quality. So only r1, r2, r3 contribute -> depth 3.
    assert pysamstats_avg == 3

    # count_coverage's own default quality_threshold=15 disagrees: it also
    # drops r3 for low base quality, leaving only r1, r2 -> depth 2.
    default_avg = _count_coverage_average(aln, chrom, start, stop)
    assert default_avg == 2
    assert default_avg != pysamstats_avg

    # Tuned: quality_threshold=0 restores r3. Default read_callback='all'
    # already matches pysamstats on duplicate/secondary/qcfail/unmapped
    # exclusion, so no custom callback is needed.
    tuned_avg = _count_coverage_average(aln, chrom, start, stop, quality_threshold=0)
    assert tuned_avg == pysamstats_avg == 3


def test_tuned_count_coverage_matches_pysamstats_across_realistic_windows(tmp_path):
    """Broader sweep: ~80x coverage over a 5kb window, many overlapping
    reads of varying start positions (no low-baseq/duplicate reads here --
    that axis is covered by the test above). Confirms the *true*
    (unrounded) means agree exactly across many distinct sub-windows, i.e.
    filter-level parity holds beyond the single hand-built region above.
    """
    import random

    random.seed(11)
    bam = tmp_path / "sweep.bam"
    contig, contig_len = "chr1", 5000
    read_len = 100
    n_reads = int(contig_len * 80 / read_len)
    reads = [
        _make_read(f"r{i}", random.randint(0, contig_len - read_len), read_len, 40, 0)
        for i in range(n_reads)
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))

    mismatches = 0
    for _ in range(50):
        w = random.choice((50, 100, 250))
        start = random.randint(0, contig_len - w - 1)
        stop = start + w

        c = pysamstats.load_coverage(
            aln, chrom=contig, start=start, end=stop, truncate=True, max_depth=300000
        )
        pysamstats_true_mean = sum(int(x) for x in c.reads_all) / len(c.reads_all)

        tuned_avg = _count_coverage_average(aln, contig, start, stop, quality_threshold=0)
        if pysamstats_true_mean != tuned_avg:
            mismatches += 1

    assert mismatches == 0


# --- 2. Numeric result: NOT reconcilable ------------------------------------


def test_production_average_depth_truncates_fractional_means(tmp_path):
    """The actual blocker. Build a region whose true mean has a fractional
    part (10 positions: 5 covered by 2 reads, 5 covered by 1 -> mean 1.5),
    with no low-baseq/duplicate reads involved (filters are not the
    variable under test here).

    ``ISClipped.average_depth`` -- the real, unmodified production method,
    calling ``statistics.mean()`` on pysamstats' int32 ``reads_all`` array
    -- returns the array-dtype-coerced (truncated) result, not 1.5. A
    count_coverage-based implementation using plain float division (either
    tuned, or the commented-out pysam-only code already sitting above
    ``average_depth`` in isclipped.py) returns the correct 1.5 and so
    necessarily disagrees with what production returns today.
    """
    bam = tmp_path / "fractional.bam"
    contig, contig_len = "chr1", 200
    # Read A spans [10, 20): depth 1 across all 10 positions.
    # Read B spans [10, 15): depth +1 across the first 5 positions only.
    reads = [
        _make_read("A", 10, 10, 40, 0),
        _make_read("B", 10, 5, 40, 0),
    ]
    _build_bam(bam, contig, contig_len, reads)
    aln = pysam.AlignmentFile(str(bam))
    chrom, start, stop = contig, 10, 20

    # True mean: (5*2 + 5*1) / 10 == 1.5
    tuned_avg = _count_coverage_average(aln, chrom, start, stop, quality_threshold=0)
    assert tuned_avg == 1.5

    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")
    production_avg = isc.average_depth(chrom, start, stop)

    # This is the hard bar this ticket could not clear: production's
    # actual output (truncated by statistics.mean's numpy int32 coercion)
    # differs from the correctly-computed count_coverage float mean.
    assert production_avg == 1  # truncated, not 1.5
    assert production_avg != tuned_avg
