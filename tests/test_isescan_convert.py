"""Converting ISEScan output into an IS table (isfinder-annotation 09).

iJump never runs ISEScan -- it is x64-only, wants a library symlink workaround
and takes ~14 minutes on this genome. The operator runs it and hands over the
tab-separated results.

The two annotations are complementary on the reference genome, which is the
argument for reading it at all: ISEScan finds three copies of an element with no
ISFinder database hit whatsoever, so an element absent from ISFinder and actively
jumping is invisible to iJump otherwise. In the other direction ISEScan needs
terminal repeats plus an ORF, so the 76 bp ``IS17`` remnant is structurally
invisible to it.
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from ijump import is_annotation, is_clustering, is_table, isescan_convert
from ijump.isclipped import ISClipped

FIXTURES = Path(__file__).parent / "fixtures"
ISESCAN_TSV = FIXTURES / "isescan" / "isescan_results.tsv"


def test_locus_columns_come_straight_from_the_isescan_row():
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    first = loci.iloc[0]
    assert first["contig"] == "NODE_1_length_3909467_cov_533.478_ID_22129"
    assert first["start"] == "708843"
    assert first["stop"] == "709881"
    assert first["family"] == "IS5"


def test_every_isescan_row_becomes_a_locus():
    """13 calls in, 13 rows out. This back-end converts; it does not filter."""
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    assert len(loci) == len(ISESCAN_TSV.read_text().strip().split("\n")) - 1


def test_names_are_the_isescan_cluster_plus_a_copy_number():
    """ISEScan reports no element name, only a family and its own cluster id.
    The id plus a copy number gives the unique per-locus name the IS table wants,
    in the ``<element>_<n>`` shape the other back-ends use -- so ``base_is_name``
    strips the counter and leaves the element.
    """
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    names = loci["is_name"].tolist()
    assert names.count("IS701_225_1") == 1
    assert names.count("IS701_225_2") == 1
    assert names.count("IS701_225_3") == 1
    assert len(set(names)) == len(names)


def test_an_element_isfinder_does_not_know_keeps_its_isescan_family():
    """The sensitivity argument for this back-end: ``new_269`` has no ISFinder
    hit at all, so ISFinder-only annotation cannot see it jumping."""
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    new = loci[loci["is_name"].str.startswith("new_269")]
    assert len(new) == 3
    assert set(new["family"]) == {"new"}


def test_the_isfinder_group_column_is_left_empty():
    """ISEScan reports a family and its own numeric cluster, not an ISFinder
    group. Putting the cluster id in the group column would file a different kind
    of value under a name that already means something.
    """
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    assert set(loci["group"]) == {""}
    assert set(loci["pident"]) == {""}


def test_the_frame_carries_every_is_table_column():
    loci = isescan_convert.read_isescan(ISESCAN_TSV)

    assert list(loci.columns) == list(is_table.COLUMNS)


def test_a_file_without_the_expected_columns_is_refused(tmp_path):
    """ISEScan's `.out` is a fixed-width report and its `.csv` a comma-separated
    twin; neither is the `.tsv` this reads. Saying so beats a KeyError.
    """
    wrong = tmp_path / "wrong.tsv"
    wrong.write_text("a\tb\tc\n1\t2\t3\n")

    with pytest.raises(isescan_convert.NotISEScanOutput) as excinfo:
        isescan_convert.read_isescan(wrong)

    assert "seqID" in str(excinfo.value)


# Clustering runs an all-vs-all blastn over the extracted loci, so these need
# BLAST+ -- the same condition the parser tier carries. Their inputs are all
# committed.
needs_blast = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ is not installed",
)

REFERENCE_GZ = FIXTURES / "isfinder" / "reference.fna.gz"


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    path = tmp_path_factory.mktemp("reference") / "reference.fna"
    with gzip.open(REFERENCE_GZ, "rb") as packed:
        path.write_bytes(packed.read())
    pysam.faidx(str(path))
    return path


@pytest.fixture(scope="module")
def converted(reference):
    return isescan_convert.convert(ISESCAN_TSV, reference)


@needs_blast
def test_isescan_calls_one_element_where_iJump_clusters_two(converted):
    """The two words "cluster" do not mean the same thing, and this back-end
    takes iJump's.

    ISEScan files all three ``new_269`` calls under one of its own clusters. Read
    off the sequences, they are not one element by iJump's rule: ``new_269_1``
    (2158 bp) aligns to ``new_269_3`` (5404 bp) over its whole length but at only
    90.7% identity, under the 95% that says a clipped read could not tell them
    apart, and it aligns to ``new_269_2`` (2315 bp) not at all. ``new_269_2`` and
    ``new_269_3`` do meet the rule, at 96.4% over 88% of the shorter.

    So the converter emits two clusters, disambiguated ``.a`` and ``.b``. That is
    the documented rule doing its job on ISEScan's calls, not a disagreement to
    paper over -- an operator who believes ISEScan can edit the column.
    """
    new = converted[converted["is_name"].str.startswith("new_269")]
    by_name = dict(zip(new["is_name"], new["cluster"]))

    assert by_name["new_269_2"] == by_name["new_269_3"]
    assert by_name["new_269_1"] != by_name["new_269_3"]
    assert set(by_name.values()) == {"new_269.a", "new_269.b"}


@needs_blast
def test_the_two_copies_of_one_element_share_a_cluster(converted):
    """``IS701_225`` at 952634 and 1719213 are two copies of one element, which
    the ISFinder back-end also puts in one cluster (as ``ISAba11``)."""
    is701 = converted[converted["is_name"].str.startswith("IS701_225")]

    assert len(set(is701["cluster"])) == 1


def test_clustering_goes_through_the_shared_core_with_the_shared_defaults(monkeypatch):
    """Same rules and thresholds as the other back-ends because it is the same
    function. Watched at the call rather than re-derived from the output, so the
    test fails if this back-end ever grows its own copy of the rules.
    """
    seen = {}

    def spy(table, reference, identity, coverage):
        seen.update(reference=reference, identity=identity, coverage=coverage)
        return table

    monkeypatch.setattr(is_annotation, "annotate_and_cluster", spy)

    isescan_convert.convert(ISESCAN_TSV, "unused_reference.fna")

    assert seen == {
        "reference": "unused_reference.fna",
        "identity": is_clustering.IDENTITY_DEFAULT,
        "coverage": is_clustering.COVERAGE_DEFAULT,
    }


@needs_blast
def test_every_locus_lands_in_a_cluster_and_carries_the_origin_flag(converted):
    assert not converted["cluster"].eq("").any()
    assert set(converted["wraps_origin"]) <= {"yes", "no"}


@needs_blast
def test_the_converted_table_is_one_ijump_run_accepts(converted, tmp_path):
    """The criterion behind "runs through both estimation modes": every check
    ``ijump run`` makes of its ``--isel`` table before touching an alignment.
    """
    path = tmp_path / "ISTable_processing.txt"
    is_table.write_is_table(converted, path)

    reread = is_table.read_is_table(path)
    assert list(reread.columns) == list(is_table.COLUMNS)
    # Both modes group by this and refuse a table without it (tickets 06, 07).
    clusters = is_table.cluster_by_name(reread)
    assert len(clusters) == len(reread)

    isc = ISClipped(FakeIsescanAlignment(), "ref", "unused.gff", str(tmp_path))
    isc.iscollect(str(path))
    isc.set_is_boundaries(100)
    assert isc.boundaries


class FakeIsescanAlignment:
    """The two contigs the committed ISEScan results call loci on."""

    references = (
        "NODE_1_length_3909467_cov_533.478_ID_22129",
        "NODE_2_length_148137_cov_371.33_ID_22131",
    )
    lengths = (3909467, 148137)


@needs_blast
def test_the_subcommand_writes_a_table_the_reader_accepts(reference, tmp_path):
    """Through the CLI and back off disk, which none of the tests above do."""
    result = subprocess.run(
        [
            "ijump",
            "isescan-convert",
            "-i",
            str(ISESCAN_TSV),
            "-r",
            str(reference),
            "-o",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    written = is_table.read_is_table(tmp_path / "ISTable_processing.txt")

    assert list(written.columns) == list(is_table.COLUMNS)
    assert len(written) == 13
    # The point of the whole back-end: a table `ijump run` will take.
    assert len(is_table.cluster_by_name(written)) == 13
