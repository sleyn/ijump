#!/usr/bin/env python3
"""Top-level ``ijump`` console-script dispatcher.

Before this module existed, iJump shipped as several independent scripts
that users invoked directly (``python3 ijump.py ...``, ``python3
combine_results.py ...``, and so on). This module wires those same
entry points up as subcommands of a single installed ``ijump`` command,
mirroring conventions like ``git`` or ``uv``:

    ijump run                 -> ijump.ijump:main()
    ijump combine-results     -> ijump.combine_results:main()
    ijump isfinder-db-parse   -> ijump.isfinder_db_parcer:main()

Each target module's ``main()`` still builds and parses its own argparse
(or getopt) command line exactly as it did as a standalone script. This
dispatcher only decides *which* module to hand control to; it does not
reinterpret or restrict that module's flags in any way. To make that work
for every target unmodified, everything typed after the subcommand name is
handed to the target as ``sys.argv[1:]`` (with ``sys.argv[0]`` set to a
friendly ``ijump <subcommand>`` label for usage/help text), then restored
once the target returns.
"""

import argparse
import importlib
import sys

# Subcommand name -> dotted module path of the script it replaces.
_SUBCOMMAND_MODULES = {
    "run": "ijump.ijump",
    "combine-results": "ijump.combine_results",
    "isfinder-db-parse": "ijump.isfinder_db_parcer",
}

# Short help text shown in `ijump --help`'s subcommand listing.
_SUBCOMMAND_HELP = {
    "run": "Run the main iJump pipeline to detect IS element rearrangements.",
    "combine-results": "Combine per-sample iJump report tables into one summary table.",
    "isfinder-db-parse": (
        "Parse a BLAST (outfmt 6) search of a genome against the ISFinder database."
    ),
}


def build_arg_parser():
    """Build a parser used only for top-level ``ijump --help``/usage output
    and for rejecting an unknown subcommand. Actual dispatch below does NOT
    route through this parser's parse_args(): argparse subparsers combined
    with a REMAINDER positional mishandle an option (e.g. ``-h``) that
    appears immediately after the subcommand name (a known argparse
    quirk -- https://bugs.python.org/issue17050), which would otherwise
    swallow each target module's own ``-h``/``--help`` instead of passing
    it through.
    """
    parser = argparse.ArgumentParser(
        prog="ijump",
        description="iJump: search for Insertion Sequence (IS) element rearrangements "
        "in evolved population sequencing data.",
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    for name in _SUBCOMMAND_MODULES:
        subparsers.add_parser(name, help=_SUBCOMMAND_HELP[name], add_help=False)

    return parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)

    # No subcommand, or a request for top-level help: let argparse print
    # the usual usage/help text (and error out on an unknown subcommand).
    if not argv or argv[0] in ("-h", "--help") or argv[0] not in _SUBCOMMAND_MODULES:
        return build_arg_parser().parse_args(argv)

    subcommand, target_args = argv[0], argv[1:]
    module = importlib.import_module(_SUBCOMMAND_MODULES[subcommand])

    # Present the target module's own argv, as if `ijump <subcommand>` were
    # itself the script being invoked directly -- so every flag/default it
    # defines (including its own -h/--help) keeps behaving exactly as it
    # did as a standalone script. Restored afterwards so callers of main()
    # (e.g. tests) aren't left with sys.argv mutated.
    original_argv = sys.argv
    sys.argv = [f"ijump {subcommand}", *target_args]
    try:
        return module.main()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    main()
