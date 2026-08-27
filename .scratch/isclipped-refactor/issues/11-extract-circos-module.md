# 11 — Extract Circos file generation into a standalone `circos` module

Status: ready-for-agent
Blocked by: —

## Why

`ISClipped.create_circos_files` (`isclipped.py:1021-1125`, ~105 lines) is an unrelated concern
riding along on the detection pipeline's seam: it consumes already-finished results
(`report_table`, `sum_by_region`, `is_coords`, `ref_len`) and writes files in a format an
external tool (Circos) reads — it has nothing to do with finding IS rearrangements. Every
`ISClipped` instance carries its colour-palette state
(`_cirocs_colors`/`_ref_colours`/`_is_colours`, `isclipped.py:187-199`) and `data_folder`
whether or not `--circos` was passed. Deletion test: delete `create_circos_files` and its
supporting state, and the detection pipeline (`ISClipped.run()`) is untouched — confirms it
was never load-bearing for detection.

Architecture review candidate 02 (`/improve-codebase-architecture`, 2026-08-11), grilled to a
shared understanding the same session. Vocabulary: `/codebase-design`.

Two things found during grilling, both load-bearing for the scope below:

1. `_rand_str` (`isclipped.py:943-945`, `@staticmethod`) and `self.session_id`
   (`isclipped.py:201`, set once in `__init__`, never read anywhere) are dead code — confirmed
   by grep, no callers.
2. `self.cutoff` (`isclipped.py:169`, default `0.005`) is read inside
   `create_circos_files` (`isclipped.py:1051`, `1065`, `1088`) to threshold which rows get
   drawn. It's a parameter this function needs, not derived pipeline state.

## Scope

### 1. Create the module

Create `circos.py` at the repo root with:

```python
def write_files(report_table, sum_by_region, is_coords, ref_len, data_folder, cutoff) -> None:
    ...
```

Move the body of `create_circos_files` (`isclipped.py:1021-1125`) into it, parameterized as
above instead of reading `self.*`. Colour assignment
(`_cirocs_colors`/`_ref_colours`/`_is_colours`) becomes local to `write_files` — nothing
outside Circos rendering reads those dicts today. Preserve the body verbatim otherwise
(relocation, not a rewrite) — including the existing `_cirocs_colors` name/typo and the
`circos.conf` template read via `script_folder + '/circos.conf'` (unaffected by the move,
`__file__`-relative path stays correct from the new module).

### 2. Delete dead code

Delete `_rand_str` (`isclipped.py:943-945`) and `self.session_id`
(`isclipped.py:201`, plus its docstring comment at `isclipped.py:200`). Neither moves to the
new module — nothing in the codebase reads either today, and carrying unused code into the
extraction just relocates it instead of removing it.

### 3. Update the call site

`ijump.py`'s `main()` (`ijump.py:143-145`) today calls `is_processing.create_circos_files()`
after `run()` returns, gated on `args.circos is True and args.estimation_mode ==
EstimationMode.AVERAGE`. Change it to call the module function directly against
`is_processing`'s attributes, same gating condition:

```python
if args.circos is True and args.estimation_mode == EstimationMode.AVERAGE:
    circos.write_files(
        is_processing.report_table, is_processing.sum_by_region,
        is_processing.is_coords, is_processing.ref_len,
        is_processing.data_folder, is_processing.cutoff,
    )
```

Add `import circos` to `ijump.py`. Remove `create_circos_files` from `ISClipped` once moved —
**no forwarding wrapper method** left behind; `ijump.py` already reaches into
`is_processing.*` attributes directly elsewhere (e.g. `result.insertions_found`), so this is
consistent with the existing call-site style.

`self.data_folder` (`isclipped.py:199`, default `'./ijump_data/'`) and `self.cutoff`
(`isclipped.py:169`) stay on `ISClipped` — they're still set/used elsewhere
(`self.cutoff` also gates report generation logic outside this ticket's scope) and are simply
passed as plain arguments at the call site.

## Out of scope

- Any change to Circos output format, column layout, or the existing `_cirocs_colors` typo.
  Byte-for-byte identical output, verified below.
- `_ref_colours`/`_is_colours` as return values or any other externally-visible artifact —
  nothing downstream reads them today; keeping them local-only is not a behavior change.

## Verification

**New characterization test** — `tests/test_circos.py`. No existing test touches Circos output
at all, so this is first coverage of the path, not a preservation of an existing test. Hand-
construct small literal `report_table`/`sum_by_region`/`is_coords`/`ref_len` fixtures (same
style as `tests/test_estimate_frequencies.py`), run today's `ISClipped.create_circos_files`
against them before the move, assert on the written file contents (`karyotype.txt`, `text.txt`,
`links.txt`, `histogram.txt`, `depth.txt`, `circos.conf`) as a golden value. After the move,
confirm the same assertions pass unchanged against `circos.write_files`.

**Manual real-sample check** — run the full pipeline with `--circos` on a real sample in
average mode before and after this change, diff the generated `ijump_data/` contents
byte-for-byte. Must match exactly.

## Done when

- `circos.py` exists with a public `write_files(...)` function; `ISClipped` no longer defines
  `create_circos_files`, `_rand_str`, `_cirocs_colors`/`_ref_colours`/`_is_colours`, or
  `session_id`.
- `ijump.py` calls `circos.write_files(...)` directly, gated on the same condition as today.
- `tests/test_circos.py` exists and passes.
- Manual real-sample before/after diff of `ijump_data/` output is clean.
- `pytest` passes from a clean clone.

## Comments
