"""Parser-level characterization golden (isfinder-annotation 01).

Pins the IS table that ``ijump isfinder-db-parse`` produces from the committed
ISFinder BLAST output, byte for byte. Both files are a few kilobytes, so this
tier runs anywhere -- no large inputs, no BLAST+ install.

This is characterization, not specification: the golden records what the parser
does today, including the behaviour the rest of this ticket series sets out to
change (the ``_\\d+`` copy suffix invented here). When a later ticket changes that
on purpose, re-pin with ``python tests/regenerate_goldens.py parser`` and review
the diff.
"""

import golden_support


def test_blast_output_produces_committed_is_table(tmp_path):
    result = golden_support.run_isfinder_db_parse(tmp_path)
    assert result.returncode == 0, result.stderr

    table_name = golden_support.ISFINDER_TABLE_NAME
    produced = tmp_path / table_name
    assert produced.exists(), f"{table_name} was not written"

    golden = golden_support.ISFINDER_GOLDEN_DIR / table_name
    assert produced.read_bytes() == golden.read_bytes(), (
        f"{table_name} differs from the golden.\n"
        f"produced:\n{produced.read_text()}\n"
        f"golden:\n{golden.read_text()}\n"
        "If the change is intended: python tests/regenerate_goldens.py parser"
    )
