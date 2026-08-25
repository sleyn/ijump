"""The IS table format: subject-id parsing, and reading both table generations.

The table ``isfinder-db-parse`` writes and ``ijump run`` reads is a headered TSV
(isfinder-annotation 03). Tables written before the header existed -- four
whitespace-separated columns, no annotation -- are still read, so an upgrade does
not strand a working setup.
"""

import pandas as pd
import pysam
import pytest

from ijump import is_table
from ijump.isclipped import ISClipped


@pytest.mark.parametrize(
    "sseqid,expected",
    [
        # The ordinary shape: name, family, group.
        ("ISAba18_IS3_IS51", ("ISAba18", "IS3", "IS51")),
        ("ISAba1_IS4_IS10", ("ISAba1", "IS4", "IS10")),
        # 11 database names carry an underscore of their own. Parsing from the
        # right keeps them whole; parsing from the left truncated them.
        ("ISBj2_B_IS5_IS5", ("ISBj2_B", "IS5", "IS5")),
        ("IS1247a_ISPpu12_IS3_IS3", ("IS1247a_ISPpu12", "IS3", "IS3")),
    ],
)
def test_subject_id_splits_from_the_right(sseqid, expected):
    assert is_table.parse_subject_id(sseqid) == expected


@pytest.mark.parametrize(
    "sseqid,expected",
    [
        ("ISAba18", ("ISAba18", "", "")),
        ("ISAba18_IS3", ("ISAba18_IS3", "", "")),
    ],
)
def test_subject_id_without_family_and_group_keeps_the_whole_id_as_name(sseqid, expected):
    """Fewer than three fields is not a name to be truncated -- it is an id with
    no annotation to recover, so nothing is guessed."""
    assert is_table.parse_subject_id(sseqid) == expected


def test_reads_a_headered_table(tmp_path):
    table = tmp_path / "ISTable_processing.txt"
    table.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\n"
        "ISAba18_1\tNODE_2\t93700\t95008\tIS3\tIS51\tISAba18\t100\n"
        "ISAba1_1\tNODE_2\t112397\t113576\tIS4\tIS10\tISAba1\t99.5\n"
    )

    read = is_table.read_is_table(table)

    assert list(read.columns) == list(is_table.COLUMNS)
    assert read.loc[0, "is_name"] == "ISAba18_1"
    assert read.loc[0, "contig"] == "NODE_2"
    assert read.loc[0, "start"] == "93700"
    assert read.loc[0, "family"] == "IS3"
    assert read.loc[1, "group"] == "IS10"
    assert read.loc[0, "cluster"] == "ISAba18"
    assert read.loc[1, "pident"] == "99.5"


def test_reads_a_legacy_headerless_table_with_the_added_columns_empty(tmp_path):
    table = tmp_path / "is_coords.txt"
    table.write_text("IS1\ttiny_contig\t900\t1000\nIS2 other_contig 10 20\n")

    read = is_table.read_is_table(table)

    assert list(read.columns) == list(is_table.COLUMNS)
    assert list(read["is_name"]) == ["IS1", "IS2"]
    assert list(read["contig"]) == ["tiny_contig", "other_contig"]
    assert list(read["stop"]) == ["1000", "20"]
    assert list(read["family"]) == ["", ""]
    assert list(read["group"]) == ["", ""]
    assert list(read["cluster"]) == ["", ""]
    assert list(read["pident"]) == ["", ""]


def test_a_legacy_table_whose_first_row_is_named_is_name_is_still_data(tmp_path):
    """The header is recognised by the whole first row, not by the first cell,
    so an IS element that happens to be called ``is_name`` is not eaten as a header."""
    table = tmp_path / "is_coords.txt"
    table.write_text("is_name\tcontig\t900\t1000\n")

    read = is_table.read_is_table(table)

    assert list(read["is_name"]) == ["is_name"]
    assert list(read["start"]) == ["900"]


def test_a_column_the_reader_does_not_know_about_survives(tmp_path):
    """A header exists so a table can carry more than its reader knows about --
    a column from a later format generation, or one an operator added."""
    table = tmp_path / "ISTable_processing.txt"
    table.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\tnotes\n"
        "ISAba18_1\tNODE_2\t93700\t95008\tIS3\tIS51\tISAba18\t100\tsplit by hand\n"
    )

    read = is_table.read_is_table(table)

    assert list(read.columns) == list(is_table.COLUMNS) + ["notes"]
    assert read.loc[0, "notes"] == "split by hand"


def test_a_headered_table_missing_an_annotation_column_reads_it_as_empty(tmp_path):
    """A hand-written table may leave out a column it has nothing to say about."""
    table = tmp_path / "ISTable_processing.txt"
    table.write_text("is_name\tcontig\tstart\tstop\n" + "IS1\ttiny_contig\t900\t1000\n")

    read = is_table.read_is_table(table)

    assert list(read.columns) == list(is_table.COLUMNS)
    assert read.loc[0, "family"] == ""


def test_written_table_reads_back_unchanged(tmp_path):
    written = pd.DataFrame(
        [
            ["ISAba18_1", "NODE_2", "93700", "95008", "IS3", "IS51", "ISAba18", "100", "no", ""],
            [
                "ISBj2_B_1",
                "NODE_2",
                "1",
                "2",
                "IS5",
                "IS5",
                "ISBj2_B",
                "98.2",
                "yes",
                "ISBj2_B_origin1",
            ],
        ],
        columns=list(is_table.COLUMNS),
    )
    table = tmp_path / "ISTable_processing.txt"
    is_table.write_is_table(written, table)

    pd.testing.assert_frame_equal(is_table.read_is_table(table), written)


def test_iscollect_keeps_coordinates_to_three_fields_for_both_formats(tmp_path, fixtures_dir):
    """``is_coords`` is a coordinate triple whatever the table's width -- Circos
    unpacks it as one, so the annotation columns must not leak into it."""
    headered = tmp_path / "ISTable_processing.txt"
    headered.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\n"
        "IS1_1\ttiny_contig\t900\t1000\tIS3\tIS51\tIS1\t100\n"
    )
    legacy = fixtures_dir / "is_coords.txt"

    for table, name in ((headered, "IS1_1"), (legacy, "IS1")):
        aln = pysam.AlignmentFile(str(fixtures_dir / "tiny.bam"))
        isc = ISClipped(aln, "unused.fna", "unused.gff", "unused_wd")
        isc.iscollect(str(table))

        assert isc.is_coords[name] == ["tiny_contig", "900", "1000"]

    assert isc.is_table is not None


def test_cluster_by_name_maps_every_row(tmp_path):
    """The cluster column is the grouping key precise mode pairs on, so it is
    read back as a plain ``is_name -> cluster`` lookup."""
    table = tmp_path / "ISTable_processing.txt"
    table.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\n"
        "IS17_1\tNODE_2\t147994\t148137\tIS5\tIS903\tISAba12\t98.6\n"
        "IS17_2\tNODE_2\t2\t77\tIS5\tIS903\tISAba12\t100\n"
        "ISAba12_1\tNODE_1\t2983841\t2984879\tIS5\tIS903\tISAba12\t98.5\n"
        "ISAba1_1\tNODE_2\t112397\t113576\tIS4\tIS10\tISAba1\t100\n"
    )

    clusters = is_table.cluster_by_name(is_table.read_is_table(table))

    assert clusters == {
        "IS17_1": "ISAba12",
        "IS17_2": "ISAba12",
        "ISAba12_1": "ISAba12",
        "ISAba1_1": "ISAba1",
    }


def test_cluster_by_name_rejects_a_legacy_table_naming_the_migration_subcommand(fixtures_dir):
    """A four-column table carries no cluster to group on. The error has to hand
    the operator the remedy, or it strands a working setup."""
    legacy = is_table.read_is_table(fixtures_dir / "is_coords.txt")

    with pytest.raises(is_table.MissingClusterColumn) as excinfo:
        is_table.cluster_by_name(legacy)

    assert is_table.MIGRATE_SUBCOMMAND in str(excinfo.value)


def test_cluster_by_name_rejects_a_row_with_the_cluster_left_blank(tmp_path):
    """An operator editing the table by hand can blank one cell. Guessing what
    that row belongs to is exactly what this ticket removes."""
    table = tmp_path / "ISTable_processing.txt"
    table.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\n"
        "IS17_1\tNODE_2\t147994\t148137\tIS5\tIS903\tISAba12\t98.6\n"
        "ISAba1_1\tNODE_2\t112397\t113576\tIS4\tIS10\t\t100\n"
    )

    with pytest.raises(is_table.MissingClusterColumn) as excinfo:
        is_table.cluster_by_name(is_table.read_is_table(table))

    assert "ISAba1_1" in str(excinfo.value)
    assert is_table.MIGRATE_SUBCOMMAND in str(excinfo.value)


def test_cluster_by_name_rejects_duplicate_names(tmp_path):
    """A dict keyed on the name would keep the last row's cluster and drop the
    other without a word -- the same reason clustering rejects duplicates."""
    table = tmp_path / "ISTable_processing.txt"
    table.write_text(
        "is_name\tcontig\tstart\tstop\tfamily\tgroup\tcluster\tpident\n"
        "IS17_1\tNODE_2\t147994\t148137\tIS5\tIS903\tISAba12\t98.6\n"
        "IS17_1\tNODE_1\t2\t77\tIS5\tIS903\tISAba53\t100\n"
    )

    with pytest.raises(ValueError, match="IS17_1"):
        is_table.cluster_by_name(is_table.read_is_table(table))
