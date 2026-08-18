# 03 — Adopt uv for dependency management, dev workflow, and PyPI publish

**What to build:** `uv` becomes the single tool for local dev (venv +
dependency sync), test running, and building/publishing the package to
PyPI. A committed `uv.lock` becomes the source of truth for exact Python
dependency versions, reused later by ticket 05's Docker image.

## Why

Grilled directly with the user: rather than uv being local-dev-only while
Docker resolves dependencies separately via pip, uv should be used
everywhere so the same locked versions apply locally, in the package build,
and inside the container — avoiding version drift between environments.

## Scope

- Depends on ticket 01's `pyproject.toml` existing.
- Run `uv lock` to generate `uv.lock` from `pyproject.toml`'s dependencies;
  commit it.
- Confirm/adjust `pyproject.toml`'s build-backend choice from ticket 01 is
  one uv's build commands (`uv build`) support (uv supports any
  PEP 517 backend — setuptools, hatchling, etc. — so this is mostly a
  compatibility check, not a required backend switch).
- Document the new dev workflow in `README.md`'s installation/development
  section: `uv sync` (create venv + install deps), `uv run pytest` (replaces
  bare `pytest` invocation), `uv run ijump ...` (replaces the ticket 02
  console script when running from a checkout rather than an installed
  package).
- Add a `uv build`/`uv publish` note to the README or a `CONTRIBUTING`-style
  doc for cutting a PyPI release (`ijump` confirmed available as a PyPI
  package name — checked via `https://pypi.org/pypi/ijump/json` returning
  404 at grilling time; re-verify at publish time in case it's been claimed
  since).
- Verify `pytest` still passes when run via `uv run pytest` (not just bare
  `pytest` in an ad-hoc venv) — this is the actual acceptance bar, since
  it's the workflow being adopted.

## Out of scope

- Actually publishing to PyPI in this ticket (this ticket makes it
  possible; a real publish is a separate, deliberate release action outside
  ticket-work scope).
- conda packaging (ticket 04) — conda's `environment.yml`/`meta.yaml` are
  independent artifacts from `uv.lock`, not derived from it (conda needs to
  express the non-Python BLAST+ dependency, which uv/PyPI cannot).
- Dockerfile changes beyond what ticket 05 scopes (ticket 05 is where the
  Docker image actually switches to `uv pip install`/`uv sync`).

## Verification

- `uv sync` succeeds from a clean clone and produces a working venv.
- `uv run pytest` passes from a clean clone.
- `uv run ijump --help` (or whichever entry point ticket 02 lands) works.
- `uv build` produces a wheel/sdist without error.
- `uv.lock` is committed and reproducible (`uv lock --check`/equivalent
  shows no drift against `pyproject.toml`).

**Blocked by:** 01 (src-layout package); soft-depends on 02 (CLI dispatch)
for the `uv run ijump` doc example, but can proceed without it if 02 hasn't
landed yet — use the pre-02 invocation form and note it needs updating.

**Status:** done, partially — blocked on `pysam`/`pysamstats` not being
resolvable/buildable via `uv` or `pip` on a current Python (see Comments;
everything else in scope is done)

- [ ] `uv.lock` generated and committed. **Not done — cannot be honestly
      generated**, see Comments. No lock file is committed.
- [x] `pyproject.toml` build-backend confirmed compatible with `uv build`.
- [x] README documents `uv sync` / `uv run pytest` / `uv run ijump` as the dev workflow (with the real limitation called out explicitly, not hidden).
- [ ] `uv run pytest` passes from a clean clone. **Not done — blocked**, see Comments.
- [x] `uv build` succeeds.

## Comments

Worked from the `refactor` branch tip (`ac2cb47`, which already has
tickets 01 and 02 merged — src-layout package and the `ijump`
console-script dispatcher). Environment: macOS/arm64, `uv 0.11.8`,
system default Python resolved by uv is 3.13.13.

**What was actually tried, in order, before concluding this is a real
blocker (not a config mistake):**

1. `uv lock` — fails immediately. Resolver is forced to `pysam==0.15.4`
   because `pysamstats==1.1.2` (the latest, and only actively-installable,
   release on PyPI — last published 2018) hard-pins `pysam<0.16` in its
   own `install_requires`. `pysam==0.15.4`'s sdist build then fails with
   `ModuleNotFoundError: No module named 'pkg_resources'` inside uv's
   isolated build env.
2. Confirmed the problem is specific to `pysamstats`' old pin, not `pysam`
   generally: `uv add pysam` (no version constraint, no `pysamstats`) in a
   throwaway scratch project resolves and installs `pysam==0.24.0`
   cleanly (it ships modern wheels).
3. Tried uv's own suggested fixes: `[tool.uv.extra-build-dependencies]
   pysam = ["setuptools<81"]` gets past the `pkg_resources` error, but
   then hits a second wall — `pysam==0.15.4`'s sdist has no bundled
   precompiled `.c` sources compatible with this build, so it needs
   `cython`, which isn't declared as a build dependency anywhere.
4. Tried `uv lock --no-build-isolation-package pysam --no-build-isolation-package pysamstats`
   — just moves the missing-build-tool problem around (now `setuptools`
   itself is missing in the non-isolated env) without fixing anything.
5. To rule out "this is a uv bug," reproduced independently with plain
   `pip` in a throwaway venv: pre-installed `setuptools<81` and `cython<3`,
   then ran `pip install --no-build-isolation pysam==0.15.4` directly.
   Still fails, now during `setup.py egg_info` itself (deeper native/
   build-configuration issue in this 2019-era sdist, unrelated to
   `pkg_resources`/Cython availability). Same result installing
   `pysamstats==1.1.2 --no-deps` under a separate Python (Homebrew
   Python 3.9, which already has a newer `pysam==0.16.0.1` installed):
   pip still tries to build `pysam==0.15.4` from source per
   `pysamstats`' own pin and hits the identical `cython`/native-build
   failure. **This confirms it's not a `uv`-specific problem** — no
   PEP 517-compliant installer can do a plain source resolve of
   `pysamstats==1.1.2` on a current Python/toolchain; it needs prebuilt
   binaries (conda) or a pinned legacy Python + build toolchain that this
   machine doesn't have.

This corroborates ticket 16's independent finding (`.scratch/isclipped-
refactor/issues/16-evaluate-dropping-pysamstats.md`, decided
`pysamstats` stays — dropping it fails on a correctness ground unrelated
to packaging) and the pre-existing README note that `pysam`/`pysamstats`
are "conda-only packages on most platforms."

**Judgment calls:**

- Did not commit a `uv.lock`. `uv lock` never produces a partial/incomplete
  lock file on failure (confirmed: no file is written at all), so there
  is nothing honest to commit. Per this ticket's own escape valve, this
  limitation is documented explicitly (README `Development setup` /
  `uv` subsection, and here) rather than worked around by, e.g., loosening
  `pyproject.toml`'s `pysam`/`pysamstats` constraints or moving them to an
  optional extra just to make `uv lock` succeed on paper — that would
  hide the real problem, not solve it, and ticket 03 is explicitly not
  the place to restructure runtime dependencies (that's ticket 16's call,
  already made, or ticket 04's conda work).
- Left `pyproject.toml`'s `[build-system]` (`setuptools.build_meta`)
  unchanged: `uv build` (see below) works with it as-is, so no backend
  switch was warranted.
- `uv build`/`uv publish` documentation was added to the README under a
  new "Releasing to PyPI" subsection rather than a separate
  `CONTRIBUTING.md`, to keep dev-workflow docs in one place; no actual
  publish was performed (out of scope) and no PyPI token was used.

**Verification results:**

- `uv sync` from a clean state: **fails** — see above. Does not produce a
  working venv.
- `uv run pytest`: **fails** — `uv run` syncs the project venv first,
  which hits the same `pysam`/`pysamstats` build failure before pytest
  ever runs. Could not distinguish pre-existing test failures from new
  ones via this path since it never gets far enough to run any tests;
  ticket 16's comments (same branch lineage, run via a conda env outside
  `uv`) record the last known-good baseline: 45 passed, 1 pre-existing
  failure (`test_read_count_mtx_rejects_invalid_orientation`, references
  a method that no longer exists, unrelated to this ticket).
- `uv run ijump --help`: **fails**, for the same reason (`uv run` syncs
  first) — even though `cli.py`'s no-subcommand path itself only imports
  `argparse`/`sys` and would work fine if the venv already existed.
- `uv build`: **succeeds** — produces `dist/ijump-1.0.4-py3-none-any.whl`
  and `dist/ijump-1.0.4.tar.gz` without touching dependency resolution.
  Verified twice (including after the README changes) to make sure
  nothing regressed it.
- `uv.lock`: not committed (see above).

**Bottom line for whoever picks up ticket 04/05:** `uv` is fully wired up
and documented as the target workflow, and works today for building the
package (`uv build`). It cannot yet manage the dev venv or run tests
end-to-end because `pysamstats==1.1.2` (PyPI's only available release)
pins an unbuildable `pysam==0.15.4`. This is the same wall ticket 04
(conda) exists to route around via prebuilt binaries — once a conda
(or otherwise binary-provided) `pysam`/`pysamstats` is available, `uv
sync --no-install-package pysam --no-install-package pysamstats` (or
similar) combined with that externally-provided pair is the likely
practical path, but that combination was not tested here since it starts
to blur into ticket 04/05's scope and this ticket's instructions were to
document the real state rather than build a workaround.

**Correction (2026-08-17, review-followups ticket 02):** The `test_read_count_mtx_rejects_invalid_orientation` failure noted above was *not* pre-existing. It was introduced by isclipped-refactor ticket 09's extraction of the read-count-matrix helper to module level in `frequency_estimation.py`, which left no delegating alias on `ISClipped`; on `master` the helper is still an `ISClipped` static method and the test passes there. It is fixed by review-followups ticket 02 (`.scratch/review-followups/issues/02-fix-read-count-mtx-test-and-baseline.md`), which repoints the test at `frequency_estimation._read_count_mtx`.
