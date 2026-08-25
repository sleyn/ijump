"""The multi-sample merge refuses samples annotated against different IS tables.

``read_reports`` joins the samples on IS identity, and an IS name is now a
cluster name -- derived from whichever IS table annotated that run. Two runs
against different tables can both produce an ``ISAba12`` that is not the same
element, and joining on the name would line them up as though it were
(isfinder-annotation 07).
"""

import pandas as pd
import pytest

from ijump import annotation_stamp
from ijump.combine_results import add_element_column, read_reports

REPORT = pd.DataFrame(
    {
        "IS Name": ["ISAba12", "ISAlw13"],
        "Annotation": ["geneA", "geneB"],
        "Chromosome": ["NODE_1", "NODE_1"],
        "Start": [100, 500],
        "Stop": [200, 600],
        "Frequency": [0.5, 0.25],
        "Depth": [100, 100],
    }
)


def _report(tmp_path, sample, fingerprint, table=REPORT):
    path = tmp_path / f"ijump_{sample}.txt"
    annotation_stamp.write_report(table, path, fingerprint)
    return str(path)


def test_samples_sharing_an_is_table_merge(tmp_path):
    """Criterion the refusal must not cost: the ordinary merge is unchanged."""
    reports = [
        _report(tmp_path, "s1", "abc123"),
        _report(tmp_path, "s2", "abc123"),
    ]

    merged = read_reports(reports, [], clonal_workflow=False, mode="average")

    assert sorted(merged.columns) == [
        "Annotation",
        "Chromosome",
        "IS Name",
        "Start",
        "Stop",
        "s1",
        "s2",
    ]
    assert len(merged) == 2
    assert merged.query('`IS Name` == "ISAba12"')["s1"].tolist() == [0.5]


def test_samples_from_different_is_tables_are_refused(tmp_path):
    reports = [
        _report(tmp_path, "s1", "abc123"),
        _report(tmp_path, "s2", "def456"),
    ]

    with pytest.raises(annotation_stamp.MixedAnnotations) as excinfo:
        read_reports(reports, [], clonal_workflow=False, mode="average")

    message = str(excinfo.value)
    assert "ijump_s1.txt" in message
    assert "ijump_s2.txt" in message


def test_the_refusal_comes_before_any_report_is_read_as_data(tmp_path):
    """A mismatch is a setup problem, so it is reported as one rather than as
    whatever the first badly-joined column happens to raise."""
    mismatched = _report(tmp_path, "s2", "def456")
    unreadable = tmp_path / "ijump_s3.txt"
    annotation_stamp.write_report(pd.DataFrame({"IS Name": ["x"]}), unreadable, "abc123")

    with pytest.raises(annotation_stamp.MixedAnnotations):
        read_reports(
            [_report(tmp_path, "s1", "abc123"), mismatched, str(unreadable)],
            [],
            clonal_workflow=False,
            mode="average",
        )


def test_the_lab_format_element_column_does_not_truncate_cluster_names():
    """The lab format collapses rows of one element together, and used to derive
    that element by stripping a numeric copy suffix off a per-locus name
    (``IS17_1`` -> ``IS17``). A cluster name has no such suffix, so the same
    regex eats a real digit: ``ISAba12``, ``ISAba11`` and ``ISAba1`` all become
    ``ISAba``, summing three elements' frequencies into one row.

    Clustering already did the collapsing, so the element is the reported name.
    """
    summary_table = pd.DataFrame({"IS Name": ["ISAba12", "ISAba11", "ISAba1", "ISAlw13", "ISVsa3"]})

    with_element = add_element_column(summary_table)

    assert with_element["Mutation"].tolist() == [
        "ISAba12",
        "ISAba11",
        "ISAba1",
        "ISAlw13",
        "ISVsa3",
    ]
