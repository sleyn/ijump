"""An ISEScan-derived IS table drives a real run, in both estimation modes.

The converter's other tests stop at the table. This is the criterion behind
"the resulting table runs through both estimation modes unmodified": hand it to
``ijump run`` on the sample alignment and see the run finish and write its
reports.

Nothing is pinned here. The ISEScan annotation calls different spans from the
ISFinder one -- deliberately, that being the argument for supporting it -- so its
outputs are not the e2e goldens and pinning them would be pinning a second
baseline nobody reads. What is asserted is that the run completes and produces
the file set for its mode.

Skips wherever the large inputs are missing, like the rest of the e2e tier.
"""

import golden_support
import pytest

from ijump import is_table, isescan_convert

pytestmark = pytest.mark.e2e

ISESCAN_TSV = golden_support.FIXTURES_DIR / "isescan" / "isescan_results.tsv"

EXPECTED_OUTPUTS = {
    "average": ("ijump_junctions.txt", "ijump_sum_by_reg.txt", "ijump_report_by_is_reg.txt"),
    "precise": ("ijump_junctions.txt", "ijump_junction_pairs.txt"),
}


@pytest.fixture(scope="module")
def isescan_is_table(tmp_path_factory):
    """The converter's output, built from the real assembly the run uses."""
    reasons = golden_support.missing_e2e_requirements()
    if reasons:
        pytest.skip("end-to-end inputs unavailable: " + "; ".join(reasons))

    reference = golden_support.e2e_data_dir() / golden_support.E2E_INPUTS["ref"]
    converted = isescan_convert.convert(ISESCAN_TSV, reference)

    path = tmp_path_factory.mktemp("isescan_table") / "ISTable_processing.txt"
    is_table.write_is_table(converted, path)
    return path


@pytest.mark.parametrize("mode", golden_support.E2E_MODES)
def test_an_isescan_table_drives_a_run(mode, isescan_is_table, tmp_path_factory):
    run_dir = tmp_path_factory.mktemp(f"isescan_e2e_{mode}")

    result = golden_support.run_e2e_pipeline(mode, run_dir, is_table_path=isescan_is_table)

    assert result.returncode == 0, f"{mode} run failed:\n{result.stdout}\n{result.stderr}"
    for filename in EXPECTED_OUTPUTS[mode]:
        assert (run_dir / "out" / filename).is_file(), f"{mode} run did not write {filename}"


def test_the_converted_table_carries_loci_the_isfinder_back_end_never_saw(isescan_is_table):
    """Why this table is worth running at all: ``new_269`` has no ISFinder hit,
    so those three loci exist in no ISFinder-derived table."""
    converted = is_table.read_is_table(isescan_is_table)
    isfinder = is_table.read_is_table(golden_support.E2E_IS_TABLE)

    assert (converted["is_name"].str.startswith("new_269")).sum() == 3
    assert not isfinder["is_name"].str.startswith("new").any()
