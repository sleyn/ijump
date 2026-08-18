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

- [ ] `ISClipped` grows a method that writes the Circos files from its own state.
- [ ] `main()` reduced to the gating condition plus one call; no attribute walk remains.
- [ ] Gating on `--circos` and average mode still lives in `main()`.
- [ ] Circos module functions unchanged — no new dependency on `ISClipped`.
- [ ] Circos output byte-identical for a `--circos` average-mode run.
- [ ] Both negative paths (no flag; precise mode) and the no-insertions path verified.

## Comments
