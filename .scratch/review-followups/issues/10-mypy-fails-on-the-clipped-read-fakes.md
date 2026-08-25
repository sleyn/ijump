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

**Status:** ready-for-agent

- [ ] `mypy src/ijump tests` reports no errors
- [ ] The fake declares the reference-position type pysam actually returns, rather than suppressing the finding
- [ ] `pre-commit run --all-files` passes with the hook installed
- [ ] `.pre-commit-config.yaml`'s comment about there being no baseline debt is true again, or rewritten to say what is actually the case
