"""Migrating a legacy IS table (isfinder-annotation 08).

An operator with a hand-curated four-column table should not have to regenerate
it because iJump started needing annotation it never carried. The migration
back-end keeps every coordinate exactly as written and fills in the rest.

Needs BLAST+, like the parser tier: it searches each locus against a database and
clusters the results. Skips cleanly without it. Its inputs are all committed --
``tests/goldens/README.md`` says what the stand-in ISFinder database is, and
``make_isfinder_db_fixture.py`` what it can and cannot stand in for.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from ijump import is_table, migrate_is_table

FIXTURES = Path(__file__).parent / "fixtures" / "isfinder"
GOLDEN = Path(__file__).parent / "goldens" / "isfinder_db_parse" / "ISTable_processing.txt"

pytestmark = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ is not installed",
)


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    """The committed masked assembly, unpacked and indexed."""
    import gzip

    import pysam

    path = tmp_path_factory.mktemp("reference") / "reference.fna"
    with gzip.open(FIXTURES / "reference.fna.gz", "rb") as packed:
        path.write_bytes(packed.read())
    pysam.faidx(str(path))
    return path


@pytest.fixture(scope="module")
def database(tmp_path_factory):
    """The stand-in ISFinder database, as a BLAST database."""
    return _blast_database(
        FIXTURES / "isfinder_db.fna", tmp_path_factory.mktemp("isfinder_db") / "isfinder"
    )


def _blast_database(fasta, prefix):
    subprocess.run(
        ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix)],
        check=True,
        capture_output=True,
    )
    return prefix


@pytest.fixture(scope="module")
def migrated(reference, database):
    legacy = is_table.read_is_table(FIXTURES / "legacy_is_table.txt")
    return migrate_is_table.migrate(legacy, reference, database)


@pytest.fixture(scope="module")
def golden():
    return is_table.read_is_table(GOLDEN)


def test_coordinates_and_names_are_carried_through_untouched(migrated):
    """The operator curated these by hand. This back-end annotates; it does not
    re-call loci, so every locus column comes out exactly as it went in."""
    legacy = is_table.read_is_table(FIXTURES / "legacy_is_table.txt")

    pd.testing.assert_frame_equal(
        migrated[["is_name", "contig", "start", "stop"]],
        legacy[["is_name", "contig", "start", "stop"]],
    )


def test_family_and_group_are_recovered_from_the_database(migrated, golden):
    """A four-column file carries no family to recover, so they come from a fresh
    search of each locus against the database.

    What this can and cannot show: the stand-in database is built from these very
    loci and labelled from this very golden, so the test proves the plumbing --
    extract, search, pick the best hit, split the subject id -- and that each
    locus picks an entry carrying its own family and group. It cannot show the
    labels are biologically right, and it cannot tell the IS5/IS903 elements
    apart from each other, since IS17, ISAba12 and ISAba53 share that pair.
    Fidelity needs the real database.
    """
    pd.testing.assert_frame_equal(migrated[["family", "group"]], golden[["family", "group"]])


def test_percent_identity_is_written_the_way_the_primary_back_end_writes_it(migrated):
    """Same identity, same text. The primary back-end reads its hits through
    pandas, so its pident reaches the table as a float and is written "100.0";
    BLAST's own output says "100.000". One table, written by two back-ends,
    should not spell one number two ways.
    """
    assert migrated["pident"].tolist() == ["100.0"] * len(migrated)


def test_clusters_match_the_primary_back_end(migrated, golden):
    """Same rules, same thresholds -- because it is the same code
    (``is_annotation.annotate_and_cluster``)."""
    assert migrated["cluster"].tolist() == golden["cluster"].tolist()


def test_origin_spanning_flags_match_the_primary_back_end(migrated, golden):
    """The two ``IS17`` fragments are one copy the assembler cut at the contig
    seam, and migration says so as the primary back-end does."""
    pd.testing.assert_frame_equal(
        migrated[["wraps_origin", "element_id"]], golden[["wraps_origin", "element_id"]]
    )


def test_a_locus_with_no_database_hit_keeps_its_coordinates(reference, monkeypatch):
    """Annotation is best-effort; the table is not. A locus that matches nothing
    stays, with its family and group empty, because the operator put it there and
    this back-end does not re-call loci.

    The search is stubbed rather than pointed at a database of unrelated
    sequence: what is under test is what migrate does with a miss, and the miss
    is more reliably produced than a database guaranteed to contain no match.
    """
    monkeypatch.setattr(
        migrate_is_table,
        "search_database",
        lambda sequences, database: dict.fromkeys(sequences, migrate_is_table.UNANNOTATED),
    )
    legacy = is_table.read_is_table(FIXTURES / "legacy_is_table.txt")

    migrated = migrate_is_table.migrate(legacy, reference, "unused")

    assert len(migrated) == len(legacy)
    assert migrated["family"].tolist() == [""] * len(legacy)
    assert migrated["pident"].tolist() == [""] * len(legacy)
    # Clustering reads the sequences, not the database, so it still groups them.
    assert migrated["cluster"].tolist() == golden_clusters()


def golden_clusters():
    return is_table.read_is_table(GOLDEN)["cluster"].tolist()


def test_best_hits_keeps_the_highest_bitscore_and_marks_a_miss(tmp_path):
    """One locus with two hits keeps the better; one with none is UNANNOTATED
    rather than absent, so the table keeps its row."""
    blast_out = tmp_path / "hits.out"
    blast_out.write_text(
        "IS17_1\tISAba12_IS5_IS903\t98.6\t120\n"
        "IS17_1\tISAba53_IS5_IS903\t91.0\t300\n"
        "ISAlw13_1\tISAlw13_IS5_IS427\t91.9\t150\n"
    )

    hits = migrate_is_table._best_hits(
        blast_out, {"IS17_1": "", "ISAlw13_1": "", "ISNothing_1": ""}
    )

    assert hits["IS17_1"] == ("IS5", "IS903", "91.0")
    assert hits["ISAlw13_1"] == ("IS5", "IS427", "91.9")
    assert hits["ISNothing_1"] == migrate_is_table.UNANNOTATED


def test_the_subcommand_writes_a_table_the_reader_accepts(reference, database, tmp_path):
    with tempfile.TemporaryDirectory() as outdir:
        result = subprocess.run(
            [
                "ijump",
                "migrate-is-table",
                "-i",
                str(FIXTURES / "legacy_is_table.txt"),
                "-r",
                str(reference),
                "-d",
                str(database),
                "-o",
                outdir,
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr

        written = is_table.read_is_table(Path(outdir) / "ISTable_processing.txt")

    assert list(written.columns) == list(is_table.COLUMNS)
    # The migrated table is one `ijump run` will accept, which is the whole point.
    assert is_table.cluster_by_name(written)
