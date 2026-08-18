"""Direct unit test for ``ISClipped.average_depth`` (ticket 12).

Unlike the other seam tests in ``test_no_results_paths.py``, this exercises
``average_depth`` against a real ``pysam.AlignmentFile`` opened on the shared
tiny fixture (``tests/fixtures/tiny.bam``) -- ``FakeAlignment`` only stubs
``.references``/``.lengths`` and can't back a real ``pysamstats.load_coverage``
call. The expected mean was cross-checked independently with
``samtools depth -a -r tiny_contig:801-900 tests/fixtures/tiny.bam``.
"""

import pysam

from ijump.isclipped import ISClipped


def test_average_depth_returns_mean_coverage_for_region(fixtures_dir):
    aln = pysam.AlignmentFile(str(fixtures_dir / "tiny.bam"))
    isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")

    assert isc.average_depth("tiny_contig", 800, 900) == 2
