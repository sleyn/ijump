"""Characterization test for ticket 11: Circos file generation.

The fixture below is hand-built (synthetic) -- there is no pre-existing captured-real-data
fixture that exercises this path (the tiny end-to-end fixture in tests/fixtures/ never
reaches ``--circos`` output in the default test run). Expected file contents were pinned by
running the pre-move ``ISClipped.create_circos_files`` body against this same fixture (before
it was relocated to ``circos.write_files``) and capturing its actual output. This is
characterization, not specification -- it documents
today's behaviour, not whether that behaviour is correct (e.g. the ``text.txt`` IS markers
intentionally write the IS start coordinate twice instead of start/stop -- an existing quirk,
preserved verbatim).

The pre-move body read two things beyond report_table/sum_by_region/is_coords/ref_len:
``self._av_depth`` (needs a real ``self.aln`` BAM handle via pysamstats) and
``self.gff.ann_pos`` (a plain nested dict). Both are now explicit parameters of
``circos.write_files``, stubbed here with a fake, deterministic ``av_depth`` lookup function
and a literal ``gff_ann_pos`` dict -- rather than wired through a real BAM/GFF, since the
behaviour under test is the file-writing logic, not depth calculation or GFF parsing.
"""

import os

import pandas as pd

import ijump.circos as circos

REF_LEN = {"contig_1": 10000, "contig_2": 5000}
IS_COORDS = {"IS1": ["contig_1", "500", "600"], "IS2": ["contig_2", "1000", "1100"]}
REPORT_TABLE = pd.DataFrame(
    {
        "IS Name": ["IS1", "IS2", "IS1"],
        "Annotation": ["geneA", "geneB<>geneC", "-"],
        "Chromosome": ["contig_1", "contig_2", "contig_1"],
        "Start": [2000, 3000, 4000],
        "Stop": [2100, 3100, 4100],
        # Third row is below cutoff (0.005) and must be excluded from text/links output.
        "Frequency": [0.5, 0.01, 0.001],
        "Depth": [40, 20, 10],
    }
)
SUM_BY_REGION = pd.DataFrame(
    {
        "ann": ["geneA", "geneB<>geneC"],
        "chrom": ["contig_1", "contig_2"],
        "start": [2000, 3000],
        "stop": [2100, 3100],
        "IS1": [30, 0],
        "IS2": [0, 5],
    }
)
GFF_ANN_POS = {
    "contig_1": {"geneA": ["geneA", "contig_1", 2000, 2100]},
    "contig_2": {
        "geneB": ["geneB", "contig_2", 3000, 3100],
        # Zero-length annotation must be skipped (ann[3] - ann[2] <= 0).
        "zero_len": ["zero_len", "contig_2", 500, 500],
    },
}
CUTOFF = 0.005

AV_DEPTH_LOOKUP = {
    ("contig_1", 2000, 2100): 50.0,
    ("contig_2", 3000, 3100): 25.0,
}


def fake_av_depth(chrom, start, stop):
    return AV_DEPTH_LOOKUP[(chrom, start, stop)]


def _read(data_folder, name):
    with open(os.path.join(data_folder, name)) as f:
        return f.read()


def _expected_files(data_folder):
    return {
        "karyotype.txt": (
            "chr - contig_1 contig_1 0 10000 green\nchr - contig_2 contig_2 0 5000 red\n"
        ),
        "text.txt": (
            "contig_1 500 500 IS1 color=vvdgreen\n"
            "contig_2 1000 1000 IS2 color=vvdred\n"
            "contig_1 2000 2000 geneA\n"
            "contig_2 3000 3000 geneB\n"
            "contig_2 3000 3000 geneC\n"
        ),
        "links.txt": (
            "contig_1 500 600 contig_1 2000 2000 color=lgreen\n"
            "contig_2 1000 1100 contig_2 3000 3000 color=lred\n"
        ),
        "histogram.txt": ("contig_1 2000 2100 30.0,0.0\ncontig_2 3000 3100 0.0,10.0\n"),
        "depth.txt": ("contig_1 2000 2100 50.0\ncontig_2 3000 3100 25.0\n"),
    }


def _assert_matches_golden(data_folder):
    expected = _expected_files(data_folder)
    for name, expected_content in expected.items():
        assert _read(data_folder, name) == expected_content, name

    conf = _read(data_folder, "circos.conf")
    assert f"karyotype = {data_folder}karyotype.txt" in conf
    assert f"{data_folder}text.txt" in conf
    assert f"{data_folder}links.txt" in conf
    assert f"{data_folder}histogram.txt" in conf
    assert f"{data_folder}depth.txt" in conf
    assert "fill_color = green, red" in conf


def test_write_files_matches_pinned_golden_output(tmp_path):
    data_folder = str(tmp_path) + "/"

    circos.write_files(
        REPORT_TABLE,
        SUM_BY_REGION,
        IS_COORDS,
        REF_LEN,
        data_folder,
        CUTOFF,
        fake_av_depth,
        GFF_ANN_POS,
    )

    _assert_matches_golden(data_folder)
