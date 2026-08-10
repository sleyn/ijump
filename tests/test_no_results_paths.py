"""Seam tests for the empty-data paths in ``isclipped.py``.

Each test drives an ``ISClipped`` instance built from ``FakeAlignment``
directly to the site that used to call ``exit()`` (see
.scratch/isclipped-refactor/issues/02-*.md), without a BAM, GFF, reference,
or BLAST+ install.
"""
import pandas as pd
import pytest

from isclipped import ISClipped, NoInsertionsFound
from fake_alignment import FakeAlignment


def test_missing_blast_output_signals_no_insertions(tmp_path):
    isc = ISClipped(FakeAlignment(), 'ref', 'unused.gff', str(tmp_path))
    with pytest.raises(NoInsertionsFound):
        isc.parseblast('missing.out', 1)


def test_all_blast_hits_below_identity_signal_no_insertions(tmp_path):
    isc = ISClipped(FakeAlignment(), 'ref', 'unused.gff', str(tmp_path))

    # blast outfmt 6 has no header row, but parseblast() reads one anyway
    # (isclipped.py:485) and overwrites the columns afterwards, so the file
    # needs one throwaway row on top of the row that matters.
    blast_out_path = tmp_path / 'low_ident.out'
    row = '1\tcontig_1\t50\t20\t0\t0\t1\t20\t10\t30\t0.001\t40\n'
    blast_out_path.write_text(row * 2)

    with pytest.raises(NoInsertionsFound):
        isc.parseblast('low_ident.out', 1)


def test_all_hits_inside_is_boundaries_signal_no_insertions():
    isc = ISClipped(FakeAlignment(), 'ref', 'unused.gff', 'unused_wd')
    isc.boundaries = [[0, 100, 'start', 'IS1', 'contig_1']]
    isc.blastout_filtered = pd.DataFrame({
        'sseqid': ['contig_1'],
        'pos_in_ref': [50],
    })

    with pytest.raises(NoInsertionsFound):
        isc.make_gene_side_regions()


def test_read_count_mtx_rejects_invalid_orientation():
    with pytest.raises(ValueError):
        ISClipped._read_count_mtx(pd.DataFrame(), 'up')
