"""The end-to-end goldens' skip path (isfinder-annotation 01).

The end-to-end tier is expected to skip almost everywhere, so the machinery
deciding that has to work on machines where the tier itself never runs -- which
is exactly where a broken skip check would go unnoticed. These tests are
committable: they need none of the large inputs.
"""

import golden_support


def test_missing_inputs_are_reported_as_skip_reasons(tmp_path, monkeypatch):
    monkeypatch.setenv(golden_support.E2E_DATA_ENV_VAR, str(tmp_path))

    reasons = golden_support.missing_e2e_requirements()

    assert reasons, "an empty data directory must yield skip reasons"
    for filename in golden_support.E2E_INPUTS.values():
        assert any(filename in reason for reason in reasons), f"{filename} not reported missing"


def test_data_dir_honours_the_environment_override(tmp_path, monkeypatch):
    monkeypatch.setenv(golden_support.E2E_DATA_ENV_VAR, str(tmp_path))
    assert golden_support.e2e_data_dir() == tmp_path

    monkeypatch.delenv(golden_support.E2E_DATA_ENV_VAR)
    assert golden_support.e2e_data_dir() == golden_support.E2E_DATA_DEFAULT


def test_populated_data_directory_reports_no_missing_inputs(tmp_path, monkeypatch):
    """The negative path above passes even for a check that always cries missing.

    So pin the other direction too: with every input present, nothing is reported
    about the inputs. Reasons about BLAST+ can still legitimately appear -- those
    depend on the machine, not on the data directory -- so only input reasons are
    asserted away.
    """
    monkeypatch.setenv(golden_support.E2E_DATA_ENV_VAR, str(tmp_path))
    for filename in golden_support.E2E_INPUTS.values():
        (tmp_path / filename).touch()
    (tmp_path / (golden_support.E2E_INPUTS["aln"] + ".bai")).touch()

    reasons = golden_support.missing_e2e_requirements()

    input_reasons = [reason for reason in reasons if str(tmp_path) in reason]
    assert input_reasons == [], input_reasons


def test_any_accepted_bam_index_spelling_satisfies_the_check(tmp_path, monkeypatch):
    """pysam reads .bai and .csi, beside the BAM or replacing its extension."""
    monkeypatch.setenv(golden_support.E2E_DATA_ENV_VAR, str(tmp_path))
    for filename in golden_support.E2E_INPUTS.values():
        (tmp_path / filename).touch()
    alignment = tmp_path / golden_support.E2E_INPUTS["aln"]

    for index_name in ("Sample.bam.bai", "Sample.bam.csi", "Sample.bai", "Sample.csi"):
        index = tmp_path / index_name
        index.touch()
        reasons = [r for r in golden_support.missing_e2e_requirements() if str(tmp_path) in r]
        assert reasons == [], f"{index_name} should satisfy the index check, got {reasons}"
        index.unlink()

    reasons = [r for r in golden_support.missing_e2e_requirements() if str(tmp_path) in r]
    assert any("index" in reason for reason in reasons), (
        f"a BAM with no index at all must be reported, got {reasons}"
    )
    assert alignment.exists()
