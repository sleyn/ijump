"""Which IS table a report was built from, recorded in the report itself.

``ijump run`` writes the report tables and ``ijump combine-results`` merges
several of them into one comparative table, joining the samples on IS identity.
That join is only sound when every sample was annotated against the same IS
table: cluster names are *derived* from the loci rather than fixed labels, so
``ISAba12`` in one run and ``ISAba12`` in another are the same element only if
both runs saw the same table. Merging across annotations would line up names
that mean different things -- silently, and worse than the per-locus names the
clusters replaced, which at least did not move (isfinder-annotation 07).

Nothing in a report's own columns answers "which table produced this", so the
writer stamps the answer on and the reader checks it. The stamp is a leading
comment line, which leaves the table itself byte-for-byte where it was.
"""

import os
from typing import Dict, Mapping, Tuple, Union

import pandas as pd

# Anything open() takes.
Path = Union[str, "os.PathLike[str]"]

# The leading comment line, minus the fingerprint.
PROVENANCE_PREFIX = "# ijump-is-table:"


class MissingProvenance(Exception):
    """A report does not say which IS table it was built from."""


class MixedAnnotations(Exception):
    """Reports built from different IS tables were about to be merged."""


def write_report(table: pd.DataFrame, path: Path, fingerprint: str) -> None:
    """Write ``table`` as a TSV stamped with the IS table it came from."""
    with open(path, "w") as report_file:
        report_file.write(f"{PROVENANCE_PREFIX} {fingerprint}\n")
        table.to_csv(report_file, sep="\t", index=False)


def read_report(path: Path) -> Tuple[pd.DataFrame, str]:
    """Read a stamped report, returning its table and its IS table fingerprint.

    A report without the stamp raises rather than being read: it was written
    before clusters existed, so its IS names are called loci rather than
    elements, and merging it with a report that names elements would join names
    that do not mean the same thing.
    """
    with open(path, "r") as report_file:
        first_line = report_file.readline()

    if not first_line.startswith(PROVENANCE_PREFIX):
        raise MissingProvenance(
            f"{os.fspath(path)} does not say which IS table it was built from. "
            "Reports written before iJump reported by cluster name one column per "
            "called locus, and those names cannot be merged with names that mean "
            "whole elements. Rerun iJump on this sample to produce a report that "
            "carries the annotation it was built against."
        )

    return (
        pd.read_csv(path, sep="\t", skiprows=1),
        first_line[len(PROVENANCE_PREFIX) :].strip(),
    )


def check_one_annotation(fingerprint_by_report: Mapping[str, str]) -> None:
    """Raise unless every report was built from the same IS table.

    Takes the whole mapping rather than a pair so the message can name every
    annotation involved: with more than two samples, "these two disagree" sends
    the reader looking in the wrong place.
    """
    by_fingerprint: Dict[str, list] = {}
    for report, fingerprint in fingerprint_by_report.items():
        by_fingerprint.setdefault(fingerprint, []).append(report)

    if len(by_fingerprint) <= 1:
        return

    groups = "\n".join(
        f"  {fingerprint}: " + ", ".join(sorted(reports))
        for fingerprint, reports in sorted(by_fingerprint.items())
    )
    raise MixedAnnotations(
        "These reports were built from different IS tables, so their IS names do "
        "not mean the same elements and merging them would line up names that "
        f"disagree:\n{groups}\n"
        "Rerun the samples against one IS table, then combine them."
    )
