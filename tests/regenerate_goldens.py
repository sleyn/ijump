#!/usr/bin/env python3
"""Re-pin the characterization goldens (isfinder-annotation 01).

    python tests/regenerate_goldens.py            # both tiers
    python tests/regenerate_goldens.py parser     # parser-level only
    python tests/regenerate_goldens.py e2e        # end-to-end only

Run this only when a behaviour change is *intended*: it overwrites the
committed expectations with whatever the working tree currently produces. The
point of the goldens is that the resulting ``git diff`` under ``tests/goldens/``
is the reviewable record of what moved -- so regenerate, then read the diff.

The end-to-end tier needs the large inputs (see ``golden_support.E2E_INPUTS``)
in ``Test/``, or in the directory named by ``IJUMP_E2E_DATA``. Each mode takes a
couple of minutes.
"""

import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import golden_support as gs  # noqa: E402


def regenerate_parser_golden():
    reasons = gs.missing_parser_requirements()
    if reasons:
        raise SystemExit("cannot regenerate the parser golden:\n  " + "\n  ".join(reasons))

    with tempfile.TemporaryDirectory() as tmp:
        result = gs.run_isfinder_db_parse(tmp)
        if result.returncode != 0:
            raise SystemExit(f"isfinder-db-parse failed:\n{result.stderr}")
        gs.ISFINDER_GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        destination = gs.ISFINDER_GOLDEN_DIR / gs.ISFINDER_TABLE_NAME
        shutil.copy2(Path(tmp) / gs.ISFINDER_TABLE_NAME, destination)
        print(f"wrote {destination.relative_to(gs.REPO_ROOT)}")


def regenerate_e2e_goldens():
    reasons = gs.missing_e2e_requirements()
    if reasons:
        raise SystemExit("cannot regenerate end-to-end goldens:\n  " + "\n  ".join(reasons))

    for mode in gs.E2E_MODES:
        with tempfile.TemporaryDirectory() as tmp:
            print(f"running {mode} mode (a couple of minutes)...")
            result = gs.run_e2e_pipeline(mode, tmp)
            if result.returncode != 0:
                raise SystemExit(f"{mode} run failed:\n{result.stdout}\n{result.stderr}")

            for golden_rel, produced_rel in gs.E2E_GOLDEN_FILES[mode].items():
                destination = gs.E2E_GOLDEN_DIR / mode / golden_rel
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(Path(tmp) / produced_rel, destination)
                print(f"wrote {destination.relative_to(gs.REPO_ROOT)}")


def main(argv):
    tiers = argv[1:] or ["parser", "e2e"]
    unknown = [tier for tier in tiers if tier not in ("parser", "e2e")]
    if unknown:
        raise SystemExit(f"unknown tier(s): {', '.join(unknown)}. Expected 'parser' and/or 'e2e'.")

    if "parser" in tiers:
        regenerate_parser_golden()
    if "e2e" in tiers:
        regenerate_e2e_goldens()


if __name__ == "__main__":
    main(sys.argv)
