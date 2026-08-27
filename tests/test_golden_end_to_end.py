r"""End-to-end characterization goldens (isfinder-annotation 01).

Pins both estimation modes' report tables produced from the sample alignment.
The alignment is 840 MB, so it is not committed and
these tests **skip** -- cleanly, not fail -- wherever it is absent, which is
every machine but the maintainer's and every CI run. The expected outputs are
committed, in ``tests/goldens/e2e/``.

Each mode's pipeline is run once per session (a couple of minutes each) and the
files are then compared byte for byte. When a change to them is intended, re-pin
with ``python tests/regenerate_goldens.py e2e`` and review the diff -- that diff
is the point of this tier. The three ``IS17_1`` / ``IS17_2`` / ``ISAba12_1``
columns in ``ijump_report_by_is_reg.txt`` collapsing into one, for instance, is
what isfinder-annotation 07 should show here.

Inputs are read from ``Test/``, or from the directory named by
``IJUMP_E2E_DATA``.
"""

import golden_support
import pytest

pytestmark = pytest.mark.e2e

# Some pinned files run to thousands of lines; a full diff would bury the signal.
MAX_DIFF_LINES = 20


@pytest.fixture(scope="session")
def e2e_run(tmp_path_factory):
    """Run one estimation mode over the large inputs, at most once per session.

    Returns a callable ``mode -> run directory``; the same directory is handed
    back on repeat calls so the per-file comparisons below share one run.
    """
    runs = {}

    def _run(mode):
        reasons = golden_support.missing_e2e_requirements()
        if reasons:
            pytest.skip("end-to-end inputs unavailable: " + "; ".join(reasons))

        if mode not in runs:
            run_dir = tmp_path_factory.mktemp(f"e2e_{mode}")
            result = golden_support.run_e2e_pipeline(mode, run_dir)
            assert result.returncode == 0, f"{mode} run failed:\n{result.stdout}\n{result.stderr}"
            runs[mode] = run_dir
        return runs[mode]

    return _run


@pytest.mark.parametrize(
    ("mode", "golden_rel"),
    [
        (mode, golden_rel)
        for mode, files in golden_support.E2E_GOLDEN_FILES.items()
        for golden_rel in files
    ],
)
def test_output_matches_golden(e2e_run, mode, golden_rel):
    run_dir = e2e_run(mode)

    produced = run_dir / golden_support.E2E_GOLDEN_FILES[mode][golden_rel]
    assert produced.exists(), f"{mode} run did not write {produced.name}"

    golden = golden_support.E2E_GOLDEN_DIR / mode / golden_rel
    assert golden.exists(), (
        f"no golden at {golden}; generate it with: python tests/regenerate_goldens.py e2e"
    )

    if produced.read_bytes() != golden.read_bytes():
        pytest.fail(
            _diff_message(mode, golden_rel, golden.read_text(), produced.read_text()),
            pytrace=False,
        )


def _diff_message(mode, golden_rel, golden_text, produced_text):
    """A failure message that shows what moved, truncated to stay readable.

    The first differing lines are what a reviewer needs.
    """
    import difflib

    diff = list(
        difflib.unified_diff(
            golden_text.splitlines(),
            produced_text.splitlines(),
            fromfile=f"golden/{mode}/{golden_rel}",
            tofile=f"produced/{mode}/{golden_rel}",
            lineterm="",
        )
    )
    shown = diff[:MAX_DIFF_LINES]
    if len(diff) > MAX_DIFF_LINES:
        shown.append(f"... ({len(diff) - MAX_DIFF_LINES} more diff lines)")
    return (
        f"{mode} mode: {golden_rel} differs from the golden.\n"
        + "\n".join(shown)
        + "\nIf the change is intended: python tests/regenerate_goldens.py e2e"
    )
