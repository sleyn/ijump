# 10 — mypy fails on the clipped-read fakes

**What to build:** `mypy src/ijump tests` comes back clean, so the pre-commit hook this repo
tells every clone to install can actually run.

Six errors, all in the fake reads that stand in for pysam in the clipped-read-search tests.
They come from one cause: a reference-position list is written as three `None`s followed by
a range of ints. mypy reads the leading literal as a list of `None` and rejects the ints
concatenated onto it — but a mixed list is exactly right here, because that is what
`get_reference_positions(full_length=True)` returns: a position for each aligned base of a
read and `None` for each clipped one. The fake has simply never said so.

**This is live debt, not a nit.** `.pre-commit-config.yaml` widened mypy from `src/ijump` to
include `tests/` on the evidence that "running `mypy src/ijump tests` surfaced zero
findings, so there was no baseline debt to justify leaving `tests/` out". That was true when
it was written. It is not now, so anyone who installs the hook as `CLAUDE.md` instructs
cannot commit until this is fixed — and the comment in that config file states as fact
something that has since stopped being true.

Fix the fake rather than the call sites, and do not silence it with `# type: ignore`: the
type is a real statement about what pysam hands back, and stating it once is what stops the
next fake drifting the same way.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] `mypy src/ijump tests` reports no errors
- [x] The fake declares the reference-position type pysam actually returns, rather than suppressing the finding
- [ ] `pre-commit run --all-files` passes with the hook installed — **not verified here**; `pre-commit` is not installed in this environment. Every command the hook runs was checked on its own and the `files:` pattern confirmed against test paths. Carried by ticket 08.
- [x] `.pre-commit-config.yaml`'s comment about there being no baseline debt is true again, or rewritten to say what is actually the case

## Comments

`mypy src/ijump tests` now reports **no issues in 55 source files**.

**The fix is a declaration, not a suppression.** `fake_clipped_read.py` gained
`ReferencePositions = List[Optional[int]]`, named and commented with what it stands for —
pysam's position per aligned base, `None` per clipped one — and `FakeRead.__init__` takes it.
The three position lists in `test_clipped_read_search.py` are now module constants annotated
with that alias.

Two changes were needed, not one, and the second is the interesting half. Annotating the
parameter alone does not help: `[None, None, None] + list(range(510, 517))` fails on the `+`
before any expected type reaches it, because mypy has already settled the left operand as
`list[None]`. Written as `[None, None, None, *range(510, 517)]` the literal is one
expression and mypy joins the element types — but *unannotated* that join degrades to
`list[Any]`, which silences the error by giving up on the type rather than by stating it.
Annotated against `ReferencePositions` it checks properly and reveals as
`list[Union[int, None]]`. Both halves are therefore load-bearing: the unpacking makes the
expression checkable, the annotation makes the check mean something.

**Verification.** `mypy src/ijump tests` clean; `ruff check` and `ruff format --check` clean;
`tests/test_clipped_read_search.py` 14/14; full suite 198 passed, 4 skipped.

`pre-commit run --all-files` could **not** be run — `pre-commit` is not installed in this
environment and the hook is not installed in this clone. What was verified instead: the
three hooks' shared `files: ^(src/ijump|tests)/` pattern matches `tests/` paths and does not
match `simulation/`, `rule-tests/` or repo-root files, and every command the hooks run
(`ruff check`, `ruff format`, `mypy src/ijump tests`) passes on its own. The hook runner
itself is unexercised here.
