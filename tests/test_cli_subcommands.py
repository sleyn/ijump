"""The set of subcommands `ijump` dispatches to (ticket 02).

The ISFinder HTML scraper was removed once its upstream web site went down, so
the dispatcher must offer exactly the subcommands listed below -- no more, so a
subcommand cannot be added without saying so here -- and treat the old
`isfinder-parse` name like any other unknown one.
"""

import pytest

from ijump import cli

# How argparse renders the full choice list. Asserting on the rendered list
# (rather than each name in turn) is what makes "exactly these" a real check:
# adding or removing one changes this string.
RENDERED_CHOICES = "{run,combine-results,isfinder-db-parse,migrate-is-table,isescan-convert}"


def test_help_lists_exactly_the_supported_subcommands(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--help"])

    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert RENDERED_CHOICES in out
    assert "isfinder-parse" not in out.replace("isfinder-db-parse", "")


def test_removed_subcommand_hits_the_unknown_subcommand_error(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["isfinder-parse", "-i", "page.html"])

    assert excinfo.value.code == 2
    assert "invalid choice: 'isfinder-parse'" in capsys.readouterr().err
