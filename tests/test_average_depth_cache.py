"""``average_depth``'s cache must not outlive its pipeline (review-followups 11).

The cache was an ``@lru_cache`` on the method, which lives on the *class*: every
``ISClipped`` that ever answered a depth query stayed reachable from it, together
with its alignment handle, junction tables and depth dictionaries. Ruff's bugbear
rule names this (``B019``); it was suppressed rather than fixed.

The caching itself earns its keep — a region's depth is asked for once per IS
entry in the per-region report, and again when Circos draws — so these tests pin
both halves: it still caches, and it lets go.
"""

import gc
import weakref

import pysam

from ijump.isclipped import ISClipped


class CountingAlignment:
    """Counts fetches. pysam's AlignmentFile.fetch is read-only, so the way to
    see whether the cache was consulted is to wrap the file, not patch it."""

    def __init__(self, alignment):
        self._alignment = alignment
        self.fetches = []
        self.references = alignment.references
        self.lengths = alignment.lengths

    def fetch(self, *args, **kwargs):
        self.fetches.append(args)
        return self._alignment.fetch(*args, **kwargs)


def _pipeline(fixtures_dir):
    aln = CountingAlignment(pysam.AlignmentFile(str(fixtures_dir / "tiny.bam")))
    return ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")


def test_a_region_is_measured_once(fixtures_dir):
    """The second answer comes from the cache, not from the alignment."""
    isc = _pipeline(fixtures_dir)

    first = isc.average_depth("tiny_contig", 800, 900)
    second = isc.average_depth("tiny_contig", 800, 900)

    assert first == second
    assert len(isc.aln.fetches) == 1


def test_a_different_region_is_measured_again(fixtures_dir):
    """Cached per region, not one answer for all of them."""
    isc = _pipeline(fixtures_dir)

    first = isc.average_depth("tiny_contig", 800, 900)
    second = isc.average_depth("tiny_contig", 0, 100)

    assert first != second


def test_a_discarded_pipeline_is_collectable_without_a_gc_pass(fixtures_dir):
    """Reference counting alone has to be enough.

    A cache that merely takes part in a reference *cycle* would pass a test that
    calls ``gc.collect()`` first; this one does not, so it fails for a cycle as
    well as for the class-level cache the ticket is about.
    """
    isc = _pipeline(fixtures_dir)
    isc.average_depth("tiny_contig", 800, 900)
    dead = weakref.ref(isc)

    was_enabled = gc.isenabled()
    gc.disable()
    try:
        del isc
        assert dead() is None, "the pipeline is still reachable after being discarded"
    finally:
        if was_enabled:
            gc.enable()


def test_two_pipelines_do_not_share_a_cache(fixtures_dir):
    """A class-level cache answers one instance's query from another's work, which
    is wrong the moment the two are opened on different alignments."""
    first = _pipeline(fixtures_dir)
    first.average_depth("tiny_contig", 800, 900)

    second = _pipeline(fixtures_dir)

    second.average_depth("tiny_contig", 800, 900)

    assert len(second.aln.fetches) == 1, "the second pipeline reused the first one's cached answer"
