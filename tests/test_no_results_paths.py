"""Seam tests for the empty-data paths in ``isclipped.py``.

Each test drives an ``ISClipped`` instance built from ``FakeAlignment``
directly to the site that used to call ``exit()`` (see
.scratch/isclipped-refactor/issues/02-*.md), without a BAM, GFF, reference,
or BLAST+ install.

The BLAST-output-emptiness checks used to be ``ISClipped.parseblast`` seam
tests; ticket 10 moved ``parseblast`` into ``clipped_read_search`` (as
``_parseblast``), so they now drive that module function directly instead
of going through an ``ISClipped`` instance. (Fuller pinned-output coverage
for the same function lives in tests/test_clipped_read_search.py.)
"""

import pandas as pd
import pytest
from fake_alignment import FakeAlignment

from ijump.clipped_read_search import Boundary, _parseblast
from ijump.frequency_estimation import _read_count_mtx
from ijump.isclipped import ISClipped, NoInsertionsFound


def test_missing_blast_output_signals_no_insertions(tmp_path):
    with pytest.raises(NoInsertionsFound):
        _parseblast(str(tmp_path / "missing.out"), 1, pd.DataFrame())


def test_all_blast_hits_below_identity_signal_no_insertions(tmp_path):
    # blast outfmt 6 has no header row, but parseblast() reads one anyway
    # (clipped_read_search.py's _parseblast) and overwrites the columns
    # afterwards, so the file needs one throwaway row on top of the row that
    # matters.
    blast_out_path = tmp_path / "low_ident.out"
    row = "1\tcontig_1\t50\t20\t0\t0\t1\t20\t10\t30\t0.001\t40\n"
    blast_out_path.write_text(row * 2)

    with pytest.raises(NoInsertionsFound):
        _parseblast(str(blast_out_path), 1, pd.DataFrame())


def test_all_hits_inside_is_boundaries_signal_no_insertions():
    isc = ISClipped(FakeAlignment(), "ref", "unused.gff", "unused_wd")
    isc.boundaries = [Boundary(0, 100, "start", "IS1", "contig_1")]
    isc.blastout_filtered = pd.DataFrame(
        {
            "sseqid": ["contig_1"],
            "pos_in_ref": [50],
        }
    )

    with pytest.raises(NoInsertionsFound):
        isc.make_gene_side_regions()


def test_read_count_mtx_rejects_invalid_orientation():
    # `_read_count_mtx` is a module-level helper in `frequency_estimation`
    # (moved there by isclipped-refactor ticket 09, formerly an `ISClipped`
    # static method). It has no public wrapper: `estimate_frequencies`, the
    # module's public entry point, always calls it with a hardcoded "left" or
    # "right" and never exposes `orientation` to its own callers, so the
    # invalid-orientation guard is unreachable from the public API. Testing
    # the private helper directly is therefore the only way to cover it.
    with pytest.raises(ValueError, match='"left" or "right"'):
        _read_count_mtx(pd.DataFrame(), "up")
