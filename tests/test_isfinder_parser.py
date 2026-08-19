"""Characterization test for ``isfinder_parser.main`` (ticket 01).

``isfinder_parser.py`` built its result tables with three self-reassigning
``DataFrame.append`` calls, which pandas removed in 2.0. This test drives
``main()`` end-to-end against a small, hand-built ISFinder BLAST HTML page
and pins the resulting ``ISTable_full.txt``/``ISTable_processing.txt``
output -- shape, column order, index, and values -- against a baseline
captured by running the pre-rewrite ``.append()``-based code under
pandas 1.3.5 (the version the project used to pin). That baseline is
reproduced verbatim here, quirks and all: the "name" column being blank and
a stray, all-blank "check" column surviving into the output are pre-existing
behaviour (appending a frame indexed by "name" into a template that treats
"name" as an ordinary column never aligns the two), not something this
rewrite introduces or fixes -- see the comments in ``isfinder_parser.py``.

The fixture covers: two IS elements on one contig (exercising per-contig
accumulation and the score-descending sort), an IS element on a second
contig (exercising accumulation *across* contigs), and BLAST hits get
filtered where the e-value is too high.
"""

import sys

from ijump import isfinder_parser


def _alignment_block(anchor_text, is_name, score, evalue, start, stop):
    return (
        f'<a href="l">{anchor_text}</a> {is_name} <span class="x">whatever</span>\n'
        f" Score = {score} bits (270), Expect = {evalue}\n"
        f" Query  {start}  ATGCATGCAT  {stop}\n\n"
    )


def _isels_ge_row(name, family, group, origin):
    return (
        f'<a href="l1">{name}</a></td><td>{family}</td><td>{group}</td>'
        f'<td><a href="l2">{origin}</a></td><td><a href="l3">extra</a></td>'
        f"<td>ignored</td></tr>"
    )


def _contig_block(name, length, rows, alignments):
    header = f" {name} some description\n\nLength={length}\n"
    table = "<table><tr>" + "".join(rows) + "</tr></table>\n"
    pre = "<pre>&gt;" + "&gt;".join(alignments) + "</pre>\n"
    return header + table + pre


def _build_fixture_html():
    contig1 = _contig_block(
        "contig1",
        1000,
        [
            _isels_ge_row("IS1", "fam1", "grp1", "org1"),
            _isels_ge_row("IS2", "fam2", "grp2", "org2"),
            _isels_ge_row("IS4", "fam4", "grp4", "org4"),
        ],
        [
            _alignment_block("IS1_anchor", "IS1", 500, "1e-40", 10, 200),
            _alignment_block("IS2_anchor", "IS2", 300, "1e-35", 300, 350),
            # Filtered out: e-value above the 1e-30 cutoff.
            _alignment_block("IS4_anchor", "IS4", 100, "1e-10", 500, 550),
        ],
    )
    contig2 = _contig_block(
        "contig2",
        500,
        [_isels_ge_row("IS3", "fam3", "grp3", "org3")],
        [_alignment_block("IS3_anchor", "IS3", 400, "1e-32", 1, 100)],
    )
    return (
        "<html><body><article>\n"
        "<b>Query=</b>"
        + contig1
        + "<b>Query=</b>"
        + contig2
        + "</article><footer>after</footer></body></html>"
    )


def test_main_produces_pinned_output_shape_and_values(tmp_path, monkeypatch):
    html_path = tmp_path / "isfinder_sample.html"
    html_path.write_text(_build_fixture_html())

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["isfinder_parser.py", "-i", str(html_path)])

    isfinder_parser.main()

    full_table = (tmp_path / "ISTable_full.txt").read_text()
    processing_table = (tmp_path / "ISTable_processing.txt").read_text()

    # Pinned against a run of the pre-rewrite `.append()`-based code under
    # pandas 1.3.5. Blank "name" column and stray blank "check" column are
    # pre-existing quirks (see module comments), preserved verbatim.
    expected_full = (
        "\tcontig\tname\tfamily\tgroup\torigin\tscore\te-value\tstart\tstop\tcheck\n"
        "IS1_1\tcontig1\t\tfam1\tgrp1\torg1\t500\t1e-40\t10\t200\t\n"
        "IS2_1\tcontig1\t\tfam2\tgrp2\torg2\t300\t1e-35\t300\t350\t\n"
        "IS3_1\tcontig2\t\tfam3\tgrp3\torg3\t400\t1e-32\t1\t100\t\n"
    )
    expected_processing = (
        "IS1_1\tcontig1\t10\t200\nIS2_1\tcontig1\t300\t350\nIS3_1\tcontig2\t1\t100\n"
    )

    assert full_table == expected_full
    assert processing_table == expected_processing
