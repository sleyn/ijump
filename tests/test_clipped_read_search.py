"""Characterization test for ticket 10: ``clipped_read_search``.

The clipped-read-search cluster (``_clboundaries``, ``collect_clipped_reads``,
``_cl_read_cov_overlap``, ``_crtable_ungapped``, ``_write_cl_fasta``,
``runblast``, ``_choosecoord``, ``parseblast``) used to live on ``ISClipped``
and could not be exercised without constructing a full instance backed by a
real or faked ``pysam.AlignmentFile`` -- the same shape of gap tickets 06/08/09
closed for ``_find_pair`` and ``assess_isel_freq``.

The fixtures below (tests/fake_clipped_read.py) are hand-built (synthetic),
not observed real data. The expected values were pinned by running today's
(pre-move) ``ISClipped`` methods directly against these same fixtures --
``_clboundaries`` as a staticmethod, ``collect_clipped_reads``/``crtable_bwds``
for the forward/backward ``_crtable_ungapped`` loop, ``_write_cl_fasta``,
``_cl_read_cov_overlap``, and ``parseblast`` against a canned BLAST output
table (same pattern as tests/test_no_results_paths.py's BLAST fixtures:
``pd.read_csv`` treats outfmt 6's first row as a header, so a throwaway row
sits on top of the row that matters). This is characterization, not
specification -- it documents today's behaviour, not whether that behaviour
is correct.

After the move, this test targets ``clipped_read_search``'s module-level
functions directly (the "pure pysam parsing" and "BLAST round-trip" seams
ticket 10 calls out) and ``clipped_read_search.search`` end-to-end (the
public seam), with the same pinned values -- proving the move didn't change
behaviour. No BLAST+ install is required: the BLAST round-trip is exercised
either directly against a canned output file, or (for `search`) through a
fake ``run_blast`` that writes a canned file instead of invoking ``blastn``.
"""

import pandas as pd
import pandas.testing as pdt
import pytest
from fake_clipped_read import FakeAlignmentFetch, FakeRead

import ijump.clipped_read_search as clipped_read_search
from ijump.clipped_read_search import NoInsertionsFound, SearchResult, search

READ_LEFT = FakeRead(
    "read_left",
    "3S7M",
    [None, None, None] + list(range(510, 517)),
    "AAA" + "CCCCCCC",
    reference_name="contig_1",
)
READ_RIGHT = FakeRead(
    "read_right",
    "7M3S",
    list(range(720, 727)) + [None, None, None],
    "GGGGGGG" + "TTT",
    reference_name="contig_2",
)
READ_UNMAPPED = FakeRead(
    "read_unmapped",
    "3S7M",
    [None, None, None] + list(range(510, 517)),
    "AAAAAAAAAA",
    is_unmapped=True,
    reference_name="contig_1",
)
READ_NO_CLIP = FakeRead(
    "read_no_clip",
    "10M",
    list(range(550, 560)),
    "NNNNNNNNNN",
    reference_name="contig_1",
)
READ_LEFT_REVERSE = FakeRead(
    "read_left_reverse",
    "2S5M",
    [None, None] + list(range(330, 335)),
    "AA" + "CCCCC",
    is_reverse=True,
    reference_name="contig_1",
)

FORWARD_BOUNDARIES = [
    [400, 650, "start", "IS1", "contig_1"],
    [650, 850, "stop", "IS1", "contig_2"],
]
BACKWARD_BOUNDARIES = [
    [400, 650, "-", "-", "contig_1"],
    [650, 850, "-", "-", "contig_2"],
]


def _fake_aln():
    return FakeAlignmentFetch(
        {
            "contig_1": [READ_LEFT, READ_UNMAPPED, READ_NO_CLIP, READ_LEFT_REVERSE],
            "contig_2": [READ_RIGHT],
        }
    )


# --- Pure pysam parsing: _clboundaries ---


def test_clboundaries_left_clip():
    assert clipped_read_search._clboundaries(READ_LEFT) == [[1, 3, "left", 510]]


def test_clboundaries_right_clip():
    assert clipped_read_search._clboundaries(READ_RIGHT) == [[8, 10, "right", 726]]


def test_clboundaries_left_clip_reverse_read():
    # is_reverse doesn't change the boundary math, only the 'reverse' column
    # downstream -- included to document that _clboundaries doesn't look at it.
    assert clipped_read_search._clboundaries(READ_LEFT_REVERSE) == [[1, 2, "left", 330]]


# --- Pure pysam parsing: _cl_read_cov_overlap ---


def test_cl_read_cov_overlap_counts_non_junction_aligned_positions():
    cov = {"contig_1": {}}
    aligned_pairs = [(i, 100 + i) for i in range(10)]

    clipped_read_search._cl_read_cov_overlap(cov, aligned_pairs, "contig_1")

    assert cov == {"contig_1": {101: 1, 102: 1, 103: 1, 104: 1, 105: 1, 106: 1, 107: 1, 108: 1}}


def test_cl_read_cov_overlap_skips_short_pair_lists():
    cov = {"contig_1": {}}

    clipped_read_search._cl_read_cov_overlap(cov, [(0, 100), (1, 101)], "contig_1")

    assert cov == {"contig_1": {}}


# --- Pure pysam parsing: _write_cl_fasta ---


def test_write_cl_fasta_writes_sequences_at_or_above_min_len(tmp_path):
    cl_table = pd.DataFrame(
        {
            "ID": [0, 1, 2],
            "sequence": ["AAA", "AA", "TTT"],
        }
    )
    fasta_path = tmp_path / "cl.fasta"

    clipped_read_search._write_cl_fasta(cl_table, str(fasta_path), 2)

    assert fasta_path.read_text() == ">0\nAAA\n\n>1\nAA\n\n>2\nTTT\n\n"


def test_write_cl_fasta_drops_sequences_below_min_len(tmp_path):
    cl_table = pd.DataFrame(
        {
            "ID": [0, 1],
            "sequence": ["AAA", "AA"],
        }
    )
    fasta_path = tmp_path / "cl.fasta"

    clipped_read_search._write_cl_fasta(cl_table, str(fasta_path), 3)

    assert fasta_path.read_text() == ">0\nAAA\n\n"


# --- BLAST round-trip: _parseblast against a canned BLAST output table ---


def _blast_out_lines(rows):
    return "\n".join(rows) + "\n"


def test_parseblast_filters_by_identity_and_assigns_coordinates(tmp_path):
    cl_table = pd.DataFrame(
        {
            "ID": [0],
            "clip_position": ["left"],
        }
    )
    cl_table.index = cl_table["ID"]

    blast_out_path = tmp_path / "cl_blast.out"
    blast_out_path.write_text(
        _blast_out_lines(
            [
                # throwaway header row -- pd.read_csv treats the first outfmt-6 row
                # as a header, so it never reaches the filtering logic.
                "0\tcontig_1\t99.0\t90\t0\t0\t1\t90\t500\t589\t1e-40\t150",
                # unique best hit for qseqid 0 -- kept
                "0\tcontig_1\t98.0\t85\t1\t0\t1\t85\t510\t594\t1e-35\t140",
                # below the 75% identity threshold -- dropped
                "999\tcontig_1\t50.0\t20\t5\t0\t1\t20\t10\t30\t1e-3\t40",
            ]
        )
    )

    result = clipped_read_search._parseblast(str(blast_out_path), 1, cl_table)

    expected = pd.DataFrame(
        {
            "qseqid": [0],
            "sseqid": ["contig_1"],
            "pident": [98.0],
            "length": [85],
            "mismatch": [1],
            "gapopen": [0],
            "qstart": [1],
            "qend": [85],
            "sstart": [510],
            "send": [594],
            "evalue": [1e-35],
            "bitscore": [140],
            "pos_in_ref": [594],
            "orientation": ["right"],
        },
        index=[0],
    )

    pdt.assert_frame_equal(result, expected)


def test_parseblast_raises_when_output_file_missing(tmp_path):
    with pytest.raises(NoInsertionsFound, match="No BLAST hits were found."):
        clipped_read_search._parseblast(str(tmp_path / "missing.out"), 1, pd.DataFrame())


def test_parseblast_raises_when_all_hits_below_identity(tmp_path):
    cl_table = pd.DataFrame({"ID": [0], "clip_position": ["left"]})
    cl_table.index = cl_table["ID"]

    blast_out_path = tmp_path / "cl_blast.out"
    # outfmt 6 has no header row, but parseblast reads one anyway, so a
    # throwaway row sits on top of the row that matters.
    row = "0\tcontig_1\t50\t20\t0\t0\t1\t20\t10\t30\t0.001\t40\n"
    blast_out_path.write_text(row * 2)

    with pytest.raises(NoInsertionsFound, match="No significant BLAST hits."):
        clipped_read_search._parseblast(str(blast_out_path), 1, cl_table)


# --- The public seam: search() end-to-end, with a fake run_blast ---


def _fake_run_blast(rows):
    """Build a run_blast stub that writes canned outfmt-6 rows to out_file
    instead of invoking blastn. Records each call it receives."""
    calls = []

    def run_blast(query_file, ref_name, out_file):
        calls.append((query_file, ref_name, out_file))
        with open(out_file, "w") as f:
            f.write(_blast_out_lines(rows))

    run_blast.calls = calls
    return run_blast


def test_search_forward_matches_pinned_golden_output(tmp_path):
    run_blast = _fake_run_blast(
        [
            "0\tcontig_1\t99.0\t90\t0\t0\t1\t90\t500\t589\t1e-40\t150",  # throwaway header row
            "0\tcontig_1\t98.0\t85\t1\t0\t1\t85\t510\t594\t1e-35\t140",
            "999\tcontig_1\t50.0\t20\t5\t0\t1\t20\t10\t30\t1e-3\t40",
        ]
    )

    result = search(
        _fake_aln(),
        FORWARD_BOUNDARIES,
        "ref",
        str(tmp_path),
        direction=1,
        run_blast=run_blast,
    )

    assert isinstance(result, SearchResult)

    expected_clipped_reads = pd.DataFrame(
        {
            "ID": [0, 1, 2],
            "IS name": ["IS1", "IS1", "IS1"],
            "IS_chrom": ["contig_1", "contig_1", "contig_2"],
            "Read name": ["read_left", "read_left_reverse", "read_right"],
            "left pos": [1, 1, 8],
            "right pos": [3, 2, 10],
            "clip_position": ["left", "left", "right"],
            "junction_in_read": [510, 330, 726],
            "reverse": [False, True, False],
            "sequence": ["AAA", "AA", "TTT"],
        },
        index=pd.Index([0, 1, 2], name="ID"),
    )
    pdt.assert_frame_equal(result.clipped_reads, expected_clipped_reads)

    expected_blast_hits = pd.DataFrame(
        {
            "qseqid": [0],
            "sseqid": ["contig_1"],
            "pident": [98.0],
            "length": [85],
            "mismatch": [1],
            "gapopen": [0],
            "qstart": [1],
            "qend": [85],
            "sstart": [510],
            "send": [594],
            "evalue": [1e-35],
            "bitscore": [140],
            "pos_in_ref": [594],
            "orientation": ["right"],
        },
        index=[0],
    )
    pdt.assert_frame_equal(result.blast_hits, expected_blast_hits)

    assert result.match_lengths == [7, 5, 7]
    assert result.read_lengths == 47
    assert result.n_reads_analyzed == 5
    assert result.cl_read_cov_overlap == {}

    assert len(run_blast.calls) == 1
    query_file, ref_name, out_file = run_blast.calls[0]
    assert ref_name == "ref"
    assert query_file == str(tmp_path / "cl.fasta")
    assert out_file == str(tmp_path / "cl_blast.out")


def test_search_backward_matches_pinned_golden_output(tmp_path):
    run_blast = _fake_run_blast(
        [
            "0\tcontig_1\t99.0\t90\t0\t0\t1\t90\t500\t589\t1e-40\t150",  # throwaway header row
            "0\tcontig_1\t97.0\t80\t2\t0\t1\t80\t520\t599\t1e-30\t130",
        ]
    )

    result = search(
        _fake_aln(),
        BACKWARD_BOUNDARIES,
        "ref",
        str(tmp_path),
        direction=0,
        run_blast=run_blast,
    )

    expected_clipped_reads = pd.DataFrame(
        {
            "ID": [0, 1, 2],
            "IS name": ["-", "-", "-"],
            "IS_chrom": ["contig_1", "contig_1", "contig_2"],
            "Read name": ["read_left", "read_left_reverse", "read_right"],
            "left pos": [1, 1, 8],
            "right pos": [3, 2, 10],
            "clip_position": ["left", "left", "right"],
            "junction_in_read": [510, 330, 726],
            "reverse": [False, True, False],
            "sequence": ["AAA", "AA", "TTT"],
        },
        index=pd.Index([0, 1, 2], name="ID"),
    )
    pdt.assert_frame_equal(result.clipped_reads, expected_clipped_reads)

    expected_blast_hits = pd.DataFrame(
        {
            "qseqid": [0],
            "sseqid": ["contig_1"],
            "pident": [97.0],
            "length": [80],
            "mismatch": [2],
            "gapopen": [0],
            "qstart": [1],
            "qend": [80],
            "sstart": [520],
            "send": [599],
            "evalue": [1e-30],
            "bitscore": [130],
            "pos_in_ref": [599],
            "orientation": ["right"],
        },
        index=[0],
    )
    pdt.assert_frame_equal(result.blast_hits, expected_blast_hits)

    assert result.match_lengths == []
    assert result.read_lengths == 0
    assert result.n_reads_analyzed == 0
    assert result.cl_read_cov_overlap == {
        "contig_1": {511: 1, 512: 1, 513: 1, 514: 1, 515: 1, 331: 1, 332: 1, 333: 1},
        "contig_2": {721: 1, 722: 1, 723: 1, 724: 1, 725: 1},
    }

    query_file, ref_name, out_file = run_blast.calls[0]
    assert query_file == str(tmp_path / "cl_bwrd.fasta")
    assert out_file == str(tmp_path / "cl_blast_bwds.out")


def test_search_raises_without_calling_run_blast_when_no_clipped_reads(tmp_path):
    run_blast = _fake_run_blast([])

    with pytest.raises(NoInsertionsFound, match="No clipped reads were found."):
        search(_fake_aln(), [], "ref", str(tmp_path), direction=1, run_blast=run_blast)

    assert run_blast.calls == []


def test_search_backward_no_clipped_reads_message(tmp_path):
    run_blast = _fake_run_blast([])

    with pytest.raises(
        NoInsertionsFound,
        match="No clipped reads were found near estimated insertion sites.",
    ):
        search(_fake_aln(), [], "ref", str(tmp_path), direction=0, run_blast=run_blast)

    assert run_blast.calls == []
