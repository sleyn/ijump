# 14 — Fix `--estimation_mode` default silently no-oping the pipeline

**What to build:** Running `ijump.py` without passing `--estimation_mode` on the
command line — the common case, since `average` is meant to be the implicit
default — must actually run the average-mode pipeline: junction, region-summary,
and report files get written, exactly as if `--estimation_mode average` had
been passed explicitly.

## Why

`ijump.py:63` declares:

```python
parser.add_argument('--estimation_mode', type=str, default=EstimationMode.AVERAGE,
                     choices=list(EstimationMode), ...)
```

`EstimationMode` (`isclipped.py:26`) is a `class EstimationMode(str, Enum)`, so
`EstimationMode.AVERAGE` is itself a `str` instance. argparse's documented
behavior is: when a default value is a `str`, it is run through `type` before
being stored on the namespace — this applies even when `--estimation_mode` is
never passed on the command line. `type=str` here means `str(EstimationMode.AVERAGE)`,
which (via `Enum.__str__`) returns the literal string `"EstimationMode.AVERAGE"`,
not the intended value `"average"`.

The effect: `args.estimation_mode` ends up as the plain string
`"EstimationMode.AVERAGE"` rather than the enum member `EstimationMode.AVERAGE`.
That string is `!=` both `EstimationMode.AVERAGE` and `EstimationMode.PRECISE`
(str-Enum equality compares against `.value`, i.e. `"average"`/`"precise"`, not
against the member's `repr`-like `str()`). `ISClipped.run()`
(`isclipped.py:546-567`) branches on `if mode == EstimationMode.AVERAGE: ... elif
mode == EstimationMode.PRECISE: ...` — with neither branch matching, the entire
body is skipped. No exception is raised, so `run()` returns
`RunResult(insertions_found=True)`, `ijump.py`'s `main()` returns normally, and
the process exits 0 having written nothing beyond `reads.txt` (the one write
that happens earlier, before the mode branch).

Found during ticket 13's manual real-sample verification
(`.scratch/isclipped-refactor/issues/13-extract-region-summary-module.md`):
invoking the CLI end-to-end with no `--estimation_mode` flag produced only
`reads.txt`/`ijump.log`, while driving `ISClipped.run(EstimationMode.AVERAGE)`
directly (bypassing argparse) against the same BAM produced the full expected
output set. Confirmed pre-existing: reproduces identically on `HEAD~1`
(ticket 12), unrelated to ticket 13's changes.

## Scope

- Make the default resolve to a real `EstimationMode.AVERAGE` enum member on
  `args.estimation_mode`, regardless of whether `--estimation_mode` is passed
  explicitly. The `type=str` conversion still needs to work correctly for
  explicit `--estimation_mode average`/`--estimation_mode precise` values typed
  by a user on the command line — this is not a request to drop the `type=`
  argument outright, only to stop it from mangling the *default*.
- Do not change `EstimationMode`'s definition, member names, or `.value`
  strings — `ticket 05` (`estimation-mode-not-validated.md`) already
  established the `str, Enum` shape and the reasoning is still sound
  (`isclipped.py:21-24`'s comment explains why: a mistyped comparison should
  raise loudly rather than silently fall through a dead branch — this bug is
  a different failure mode, the *default* itself resolving to the wrong
  value, not a mistyped comparison).

## Verification

- Regression test: invoke the CLI argument parser (or `main()`) with no
  `--estimation_mode` flag and assert the resulting `args.estimation_mode`
  compares equal to `EstimationMode.AVERAGE` (and is usable as such by
  `ISClipped.run()`'s `mode ==` branches).
- Confirm `--estimation_mode average` and `--estimation_mode precise` passed
  explicitly still parse to the correct enum members (no regression on the
  explicit path).
- Manual real-sample run (reuse the ticket 13 verification command) with no
  `--estimation_mode` flag: confirm `ijump_junctions.txt`,
  `ijump_sum_by_reg.txt`, and `ijump_report_by_is_reg.txt` are now written
  (currently absent).

**Blocked by:** None — can start immediately

**Status:** ready-for-agent

- [ ] `args.estimation_mode` resolves to a genuine `EstimationMode.AVERAGE` member when `--estimation_mode` is omitted from the CLI.
- [ ] `--estimation_mode average` / `--estimation_mode precise` passed explicitly still resolve correctly.
- [ ] A regression test covers the omitted-flag case.
- [ ] Manual real-sample run with no `--estimation_mode` flag produces the full average-mode output set (junctions, region summary, report), not just `reads.txt`.
- [ ] `pytest` passes from a clean clone.

## Comments
