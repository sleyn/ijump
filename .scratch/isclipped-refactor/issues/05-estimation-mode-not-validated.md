# 05 — Invalid --estimation_mode silently matches neither branch

Status: ready-for-agent
Blocked by: —

## The bug

`ijump.py`:

```python
parser.add_argument('--estimation_mode', type=str, default=EstimationMode.AVERAGE,
                    help="Specifies how the IS frequency will be esimated. 'average' - by averaging the region coverage"
                         " and number of clipped reads. Or 'precise' - iJump will try to separate each insertion event.")
```

`--estimation_mode` accepts any string. `args.estimation_mode` is later compared against
`EstimationMode.AVERAGE` / `EstimationMode.PRECISE` (ticket 03) at every branch point —
`write_empty_outputs`, `check_junctions_presence`, the two workflow branches in `main`, and
the `circos` guard. If the value passed on the command line is neither, every one of those
comparisons is false: the run doesn't error, it just silently skips both workflow branches
and produces no meaningful output.

This is the same silent-fallthrough failure mode ticket 03 closed for internal comparisons
(a typo'd literal never matching) — but that fix only guards code that already writes
`EstimationMode.PRECISE`/`AVERAGE` literals. It does nothing for a user-supplied value, since
`argparse` stores the raw string and never checks it against the enum.

## Scope

Add `choices=list(EstimationMode)` (or equivalent) to the `--estimation_mode` argument so
`argparse` rejects an invalid value at parse time with a clear error, instead of the current
silent no-op run.

## Verification

A CLI-level test (subprocess, following the pattern in `tests/conftest.py`'s `run_ijump`
fixture) asserting a bogus `--estimation_mode` value exits non-zero with an argparse usage
error, rather than exiting 0 having done nothing.

## Comments

Filed while reviewing ticket 03 (`fix: Compare estimation mode against 'precise', not the
typo 'presice'`, commit 9aa9f6c) — flagged as out of scope there since that ticket was
specifically about the internal `'presice'`/`'precise'` comparison, not CLI input validation.

Implemented in 1edec5d (`fix: Reject invalid --estimation_mode at parse time`): added
`choices=list(EstimationMode)` to the argparse argument and a CLI-level test
(`tests/test_estimation_mode_validation.py`) asserting a bogus value exits non-zero with an
argparse "invalid choice" error.
