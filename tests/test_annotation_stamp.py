"""Reports say which IS table they were built from.

Cluster names are derived from the loci rather than fixed labels, so the same
name can mean different elements in two runs annotated against different tables
or different references. The multi-sample merge joins samples on those names, so
merging across annotations would silently misalign them -- strictly worse than
the per-locus names it replaces, which at least were stable. Each report carries
a digest of its IS table so the merge can check first
(isfinder-annotation 07).
"""

import pandas as pd
import pandas.testing as pdt
import pytest

from ijump import annotation_stamp

TABLE = pd.DataFrame(
    {
        "IS Name": ["ISAba12", "ISAlw13"],
        "Chromosome": ["NODE_1", "NODE_2"],
        "Frequency": [0.5, 0.25],
    }
)


def test_a_stamped_report_round_trips(tmp_path):
    report = tmp_path / "ijump_s1.txt"
    annotation_stamp.write_report(TABLE, report, "abc123")

    read, fingerprint = annotation_stamp.read_report(report)

    assert fingerprint == "abc123"
    pdt.assert_frame_equal(read, TABLE)


def test_the_stamp_is_the_first_line_and_a_comment(tmp_path):
    """A leading comment line keeps the table itself where it has always been, so
    anything reading the file with a `#`-aware reader is unaffected."""
    report = tmp_path / "ijump_s1.txt"
    annotation_stamp.write_report(TABLE, report, "abc123")

    lines = report.read_text().splitlines()
    assert lines[0] == "# ijump-is-table: abc123"
    assert lines[1].split("\t")[0] == "IS Name"


def test_an_unstamped_report_is_refused(tmp_path):
    """A report without the stamp predates the cluster columns, so its names mean
    called loci rather than elements and it cannot be merged with one that does.
    The message has to say what to do, not just what is missing."""
    report = tmp_path / "ijump_s1.txt"
    report.write_text("IS Name\tChromosome\tFrequency\nISAba12\tNODE_1\t0.5\n")

    with pytest.raises(annotation_stamp.MissingStamp) as excinfo:
        annotation_stamp.read_report(report)

    assert "ijump_s1.txt" in str(excinfo.value)
    assert "rerun" in str(excinfo.value).lower()


def test_reports_from_one_table_agree(tmp_path):
    first, second = tmp_path / "ijump_a.txt", tmp_path / "ijump_b.txt"
    annotation_stamp.write_report(TABLE, first, "abc123")
    annotation_stamp.write_report(TABLE, second, "abc123")

    # Does not raise.
    annotation_stamp.check_one_annotation({str(first): "abc123", str(second): "abc123"})


def test_reports_from_different_tables_are_refused(tmp_path):
    """The failure this exists to prevent: two samples whose `ISAba12` are not
    the same element, joined on the name as though they were."""
    with pytest.raises(annotation_stamp.MixedAnnotations) as excinfo:
        annotation_stamp.check_one_annotation({"ijump_a.txt": "abc123", "ijump_b.txt": "def456"})

    message = str(excinfo.value)
    assert "ijump_a.txt" in message
    assert "ijump_b.txt" in message
    assert "abc123" in message and "def456" in message
