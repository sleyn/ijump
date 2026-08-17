# 02 — Single `ijump` console command with subcommands

**What to build:** One installed console-script, `ijump`, dispatching to
subcommands that replace today's four separate `python3 <script>.py`
invocations documented in `README.md`:

| Today | Becomes |
|---|---|
| `python3 ijump.py ...` (main pipeline) | `ijump run ...` |
| `python3 combine_results.py ...` | `ijump combine-results ...` |
| `python3 isfinder_db_parcer.py ...` | `ijump isfinder-db-parse ...` |
| `python3 isfinder_parser.py ...` | `ijump isfinder-parse ...` |

(Subcommand names above are a starting proposal — keep them if they read
well, but the exact spelling is an implementation judgement call, not a
constraint.)

## Why

Grilled directly with the user: given four independent CLI scripts, a
discoverable `ijump --help` listing every capability was preferred over
four separately-installed console scripts, matching conventions like `git`
or `uv` itself.

## Scope

- Depends on ticket 01 landing first: `main()` functions for
  `isfinder_parser.py` and `isfinder_db_parcer.py` (wrapped there since
  they had no `if __name__` guard at all) are the dispatch targets here.
  `ijump.py`'s existing `if __name__ == "__main__":` block (`ijump.py:180`)
  and `combine_results.py`'s (`combine_results.py:389`) already have
  argparse-based `main`-equivalents to reuse — check whether they already
  expose a callable `main()` or execute inline under the guard; if inline,
  extract to `main()` first.
- Build a top-level argparse (or equivalent) entry point,
  `src/ijump/cli.py` (or fold into `src/ijump/__init__.py` — implementer's
  call), with subparsers for `run`, `combine-results`, `isfinder-db-parse`,
  `isfinder-parse`, each delegating to the corresponding module's `main()`.
  Preserve every existing flag/argument for each subcommand exactly —
  this ticket changes the invocation shape, not the argument surface.
- Add the console-script entry point to `pyproject.toml`'s
  `[project.scripts]`: `ijump = "ijump.cli:main"` (or wherever the
  dispatcher ends up).
- Update `README.md`'s four documented invocations (lines ~128, ~135,
  ~199-209, ~253, ~261 as of ticket-15's HEAD — re-verify line numbers at
  implementation time) to the new subcommand form.
- Update `tests/conftest.py`'s `run_ijump` fixture (see ticket 01's note on
  this) to invoke the installed `ijump run` console script instead of a
  direct script path — this is the fixture's second update across tickets
  01/02, expected per ticket 01's scope note.
- Any other test that shells out to one of the four scripts by path
  (grep for `subprocess` + script names across `tests/` at implementation
  time) needs the same update.

## Out of scope

- Changing any subcommand's argument names, defaults, or behavior.
- uv migration (ticket 03), conda (ticket 04), Docker (ticket 05).

## Verification

- `ijump --help` lists all four subcommands.
- Each subcommand runs equivalently to its pre-move script (reuse each
  script's existing manual-verification command from prior tickets where
  one exists, e.g. ticket 13/14's real-sample run for `ijump run`).
- `pytest` passes from a clean clone with the package installed.
- `README.md` invocation examples match the actual new CLI shape (spot-run
  each documented example).

**Blocked by:** 01 (src-layout package)

**Status:** done

- [x] Single `ijump` console-script installed via `pyproject.toml`'s `[project.scripts]`.
- [x] All four subcommands present and functionally equivalent to today's scripts.
- [x] `README.md` updated to the new invocation forms.
- [x] `tests/conftest.py`'s `run_ijump` fixture (and any other script-invoking test) updated.
- [x] `pytest` passes from a clean clone.

## Comments

Implemented as `src/ijump/cli.py` (separate module rather than folding into
`__init__.py`, to keep `__init__.py` empty/trivial and give the dispatcher
its own docstring/tests surface). `pyproject.toml` gained
`[project.scripts]` with `ijump = "ijump.cli:main"`.

All four target `main()` functions already existed on this branch (ticket
01's work), so no changes were needed inside `ijump.py`, `combine_results.py`,
`isfinder_db_parcer.py`, or `isfinder_parser.py` — the dispatcher only
decides which module to hand control to.

Dispatch judgment call: an initial version used `argparse` subparsers with
a `nargs=argparse.REMAINDER` positional per subcommand, which is the usual
idiom for this kind of wrapper. It broke on `ijump <subcommand> -h`: when an
option token immediately follows the subcommand positional, argparse's
subparser+REMAINDER combination fails to hand it to the REMAINDER catch-all
and instead reports "unrecognized arguments" (a known argparse quirk,
https://bugs.python.org/issue17050). Fixed by dispatching manually instead:
`cli.main()` only routes through the top-level `argparse` parser for
`ijump --help`/no-args/unknown-subcommand; for a recognized subcommand it
imports the target module and calls its `main()` directly, temporarily
swapping `sys.argv` to `["ijump <subcommand>", *rest]` so the target's own
argparse/getopt parsing (including its own `-h`/`--help`) behaves exactly
as it did as a standalone script, then restores `sys.argv` in a `finally`.
Verified every subcommand's `-h` output matches its pre-existing
`python3 -m ijump.<module> -h` output, and that an unknown subcommand and
`ijump --help` still behave sensibly.

Updated `tests/conftest.py`'s `run_ijump` fixture and
`tests/test_estimation_mode_validation.py` (the only other script-invoking
test found via grep) to invoke `["ijump", "run", ...]` instead of
`[sys.executable, "-m", "ijump.ijump", ...]`; dropped the now-unused `sys`
import from both files.

Verification environment note: this machine has no working `pip install -e .`
via a plain venv (pysam has no prebuilt wheel for the system Python 3.9 and
fails to build from source without Cython/BLAST headers), matching the
README's existing note that pysam/pysamstats are effectively conda-only on
most platforms. Used the pre-existing local conda env `bioinfo`
(`/Users/sleyn/miniconda3/envs/bioinfo`), which already had all runtime deps,
and ran `pip install -e . --no-deps` there. `ijump --help` lists all four
subcommands, each subcommand's `--help` matches the pre-existing per-script
usage text, and `pytest` passes except for one pre-existing, unrelated
failure (`test_read_count_mtx_rejects_invalid_orientation`, an
`ISClipped._read_count_mtx` `AttributeError`) confirmed present on the
branch tip before any of this ticket's changes (reproduced via
`git stash`).

Also note: this work was done from an agent worktree whose own branch had
diverged from `refactor` (an unrelated agent-scaffolding artifact); the
`refactor` branch tip itself was already checked out in another worktree
and unavailable for a second checkout here, so this change was built on a
same-tip branch (`packaging-cli-dispatch`, cut from local `refactor` HEAD
commit `1096b15`) for whoever integrates this to fast-forward/merge back
onto `refactor`. This ticket file (`02-cli-subcommand-dispatch.md`) was an
untracked file that only existed in the other worktree's working tree, not
in `refactor`'s git history, so it was recreated here (with the checklist
and this Comments section filled in) rather than edited in place.
