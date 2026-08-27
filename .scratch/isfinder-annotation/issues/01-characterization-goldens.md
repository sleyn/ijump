# 01 — Characterization goldens for the current IS-annotation path

**What to build:** A test suite that pins today's IS-annotation behaviour exactly as it
is, so that every later ticket in this series lands as a reviewable diff instead of an
assertion. Running the suite on unchanged code passes; running it after a behaviour
change fails loudly and shows what moved.

Two tiers, because the alignment used for end-to-end runs is far too large for CI:

- **Parser-level, committable:** the ISFinder BLAST output already in the repo, through
  `isfinder-db-parse`, to the IS table. Small inputs and outputs, runs anywhere.
- **End-to-end, local-only:** both estimation modes' report tables and the Circos inputs,
  generated from the sample alignment. Marked so it is skipped when the large inputs are
  absent.

The baseline was confirmed reproducible on 2026-08-24: re-running the BLAST search against
a fresh ISFinder database clone and re-running the parser reproduces the committed IS
table byte-for-byte.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] Parser-level golden: committed BLAST output produces the committed IS table, asserted byte-for-byte
- [x] End-to-end goldens: both estimation modes' report tables and the Circos input files pinned
- [x] End-to-end goldens skip cleanly (not fail, not error) when the large alignment inputs are unavailable
- [x] Goldens are regenerable by a documented command, so an intended change can be re-pinned deliberately
- [x] The whole suite passes against unmodified `master`

## Outcome (2026-08-24)

Landed as `tests/goldens/` (committed expectations, README on the tiers and the re-pin
command), `tests/golden_support.py` (input lookup, skip check, pipeline invocation — shared
by the tests and the regenerator), `tests/regenerate_goldens.py`, and three test modules:
`test_golden_isfinder_db_parse.py` (parser tier), `test_golden_end_to_end.py` (`e2e` marker),
`test_golden_skip_conditions.py` (the skip check itself, committable so it is exercised on
machines where the tier never runs).

Determinism was checked for both modes by diffing two independent runs in separate
directories: all eight average files and all three precise files byte-identical.

The two tiers are chained — the end-to-end run's `--isel` points at the parser tier's own
golden, not at a copy in the data directory, so re-pinning the parser golden feeds the new
IS table into the end-to-end tier instead of leaving it on a stale one.

**One production change was required, and it is not characterization.** Precise mode does not
run at all on the current branch: `junction_pairing.find_pairs` sorts its clusters by writing
back into its input arrays, and precise mode feeds it `Series.to_numpy()` results, which
pandas 3 hands out non-writeable. The run died with `ValueError: assignment destination is
read-only` before writing a single file — same family as the pandas-3 fallout already fixed in
`f84850d`. `find_pairs` now copies its inputs at entry, covered by two new tests in
`tests/test_find_pair.py` and committed separately so the goldens commit stays additive. So the acceptance criterion "passes against unmodified `master`"
holds for the parser tier and the average end-to-end tier, but *cannot* hold for the precise
tier: on unmodified `master` there is no precise-mode output to pin.

**Not pinned, deliberately:** `reads.txt` (a ~900 KB dump of every clipped read, not a report
table) and `ijump_data/circos.conf` (embeds the run directory's absolute path, so it differs
on every run by construction).

`tests/goldens/.gitignore` re-includes the directory: the root `.gitignore` lists bare output
filenames (`ijump_junctions.txt`, `reads.txt`, ...) which match at any depth and would
otherwise swallow the goldens.

## Review follow-ups applied

From `/code-review` (standards + spec axes) on the first commit:

- End-to-end `--isel` re-pointed at the parser golden (was `Test/ISTable_processing.txt`,
  untracked and hand-placed, which would have gone stale from ticket 03 onward).
- Precise-mode determinism checked, having only been checked for average.
- BAM index check widened from `.bam.bai` alone to any spelling pysam accepts, with a test.
- `test_golden_skip_conditions.py` gained the positive direction: a populated directory
  reports no missing inputs. Without it, a check that always cried "missing" would pass.
- `run_isfinder_db_parse` hoisted into `golden_support` — the parser invocation had been
  spelled out twice, contradicting that module's own docstring.
- `blast.out` provenance recorded in the goldens README; `conftest.py`'s duplicated
  `REPO_ROOT`/`FIXTURES_DIR` now import from `golden_support`; `_diff_message`'s unused
  `context` parameter inlined; `gs` alias spelled out.

Not applied: the standards axis flagged `**Status:** done` as outside
`docs/agents/triage-labels.md`'s five strings. It matches existing practice in this tracker
(`isclipped-refactor/16`), so the vocabulary gap is the doc's, not this ticket's.
