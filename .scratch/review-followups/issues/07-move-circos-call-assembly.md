# 07 — Stop `main()` reaching into `ISClipped` to assemble the Circos call

**What to build:** Asking for a Circos diagram is one call from the CLI layer.
Today `main()` reads eight attributes off the pipeline object — including a
two-hop walk into its GFF object — to build the argument list itself, so the CLI
entry point knows the internal layout of a class it should only be driving.

## Why

The final block of `ijump.py`'s `main()` gates on the `--circos` flag and the
average-mode check, then passes eight separate attributes of the pipeline object
into the Circos writer, one of them via `is_processing.gff.ann_pos`.

Two problems the review named. It is **Feature Envy**: the knowledge of which
pieces of pipeline state a Circos diagram needs lives in the CLI rather than on
the object holding that state. And it is a **Message Chain**: `main()` now
depends on `ISClipped` having a `gff` attribute which in turn has `ann_pos`, so
an internal change to that structure breaks the entry point.

Concretely, this is why ticket 06 exists — a parameter rename inside the Circos
module had to be chased into a call site in `main()`, and any future change to
the writer's signature will have to be too.

## Scope

- Add a method on `ISClipped` that writes the Circos files, taking whatever the
  caller genuinely knows (which is: nothing — the flag check stays outside) and
  sourcing the rest from `self`.
- Reduce `main()` to the gating condition plus that one call. **Keep the gating
  in `main()`**: whether the user asked for Circos, and whether the run is in
  average mode, are CLI-level concerns and belong where the parsed arguments
  live. This ticket moves the *assembly*, not the *decision*.
- Leave the Circos module's own functions as module-level functions with explicit
  parameters. The point is to stop the CLI from knowing the argument list, not to
  give the Circos writer a dependency on `ISClipped`. The new method is the
  adapter between them.

## Out of scope

- Reducing the Circos writer's parameter count, or introducing a type to bundle
  its arguments. That is ticket 09's territory, and the two tickets deliberately
  meet at this signature — see the note below.
- Changing when Circos files are written, or their contents.

## Verification

- A run with `--circos` in average mode produces byte-identical Circos output to
  before.
- A run without `--circos`, and a run with `--circos` in precise mode, both still
  skip Circos generation. The current gating is easy to get subtly wrong when
  moving it; check both negative paths, not just the positive one.
- The no-insertions-found early-exit path still behaves as before — it returns
  before this block today, so confirm the new method is not reachable when no
  insertions were found.
- `grep -n "is_processing\." src/ijump/ijump.py` shows the walk is gone.
- `ruff`/`mypy` clean on `src/ijump/`.

**Blocked by:** 06 — it renames one of the parameters being passed at this exact
call site, so doing this first would force ticket 06 to re-chase a site that had
just moved. Ticket 06 is small; take it first.

**Status:** ready-for-agent

- [x] `ISClipped` grows a method that writes the Circos files from its own state.
- [x] `main()` reduced to the gating condition plus one call; no attribute walk remains.
- [x] Gating on `--circos` and average mode still lives in `main()`.
- [x] Circos module functions unchanged — no new dependency on `ISClipped`.
- [x] Circos output byte-identical for a `--circos` average-mode run.
- [x] Both negative paths (no flag; precise mode) and the no-insertions path verified.

## Comments

Implemented in commit `4cf46c3` on `ticket-07-circos-call-assembly`.

- Added `ISClipped.write_circos_files(self)` in `src/ijump/isclipped.py`, calling
  `circos.write_files(...)` with the same eight positional arguments, in the same
  order, sourced from `self` (`self.report_table`, `self.sum_by_region`,
  `self.is_coords`, `self.ref_len`, `self.data_folder`, `self.cutoff`,
  `self.average_depth`, `self.gff.ann_pos`). Added `circos` to `isclipped.py`'s
  existing `from ijump import ...` line; no circular-import risk since
  `circos.py` imports only stdlib.
- `main()` in `src/ijump/ijump.py` now reads:
  ```python
  if args.circos is True and args.estimation_mode == EstimationMode.AVERAGE:
      is_processing.write_circos_files()
  ```
  The gating condition is untouched/unmoved. Removed the now-unused
  `from ijump import circos` import from `ijump.py`.
- `circos.py` (including `write_files`'s 8-parameter signature) was not touched,
  per the ticket's explicit "leave this signature exactly as-is" instruction —
  ticket 09's territory.

Verification performed:
- `ruff check src/ijump/` and `mypy src/ijump/`: both clean.
- `grep -n "is_processing\." src/ijump/ijump.py`: only remaining hits are
  `is_processing.iscollect(...)`, `is_processing.set_is_boundaries(...)`,
  `is_processing.run(...)`, and `is_processing.write_circos_files()` — the
  eight-attribute walk (including `is_processing.gff.ann_pos`) is gone.
- Full test suite (`pytest tests/`, 47 tests) passes in an env with the package
  installed as an editable console-script (`ijump-verify` conda env, since the
  repo's own environment isn't set up for `pip install -e .` in this worktree).
  Env was restored to its original (ijump-not-importable) state afterward so as
  not to disturb other concurrent worktree sessions sharing it.
- Positive path (byte-identical Circos output): `test_circos.py`'s existing
  golden-output test against `circos.write_files` is untouched and still passes
  (that function's code wasn't modified). Additionally verified directly that
  `ISClipped.write_circos_files` forwards the exact same 8 values in the exact
  same order as the old inline call (monkeypatched `circos.write_files` and
  asserted the captured `args` tuple), so the positive path is byte-identical
  by construction.
- Negative paths, driven end-to-end via the CLI against the `tests/fixtures`
  tiny fixture set (which yields no insertions found): confirmed no
  `ijump_data/` (Circos output) directory is created for (a) average mode
  without `--circos`, (b) precise mode with `--circos`, and (c) average mode
  with `--circos` but no insertions found (the early-exit-before-this-block
  path). All three returned exit code 0 with no Circos output, matching
  pre-change behavior.
- `code-review` skill run (Standards + Spec sub-agents in parallel against
  `refactor..HEAD`): both reports came back clean — no hard standards
  violations (repo has no CODING_STANDARDS.md/CONTRIBUTING.md), no missing or
  partial spec requirements, no scope creep. Standards review flagged two
  minor non-blocking judgement calls: (1) `ISClipped`'s import list and
  responsibility surface is growing (now also owns file-writing delegation)
  — worth watching across future tickets, not actionable here; (2) the
  3-line comment above `write_circos_files` is slightly more
  self-justifying than strictly necessary — left as-is since it documents
  the adapter's intent per the ticket's own framing.
