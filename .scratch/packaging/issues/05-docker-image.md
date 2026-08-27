# 05 — Dockerfile: local-source build, apt BLAST+, uv-installed package

**What to build:** A `Dockerfile` (+ `.dockerignore`) that builds a
self-contained image from the repo checkout — no dependency on a prior
PyPI publish — with BLAST+ installed via apt and the `ijump` package
installed via `uv`.

**Contingent dependency:** `.scratch/isclipped-refactor/issues/16-evaluate-dropping-pysamstats.md`
is evaluating whether `pysamstats` (a compiled, conda-friendlier-than-pip
dependency) can be dropped entirely. Check that ticket's outcome before
assuming `uv sync`/`uv pip install .` alone covers every runtime
dependency in the image — if 16 does *not* land as a swap, `pysamstats`
being conda-only per `README.md` may need its own image-install step
beyond plain `uv sync` (worth a fresh look at implementation time).

## Why

Grilled directly with the user across two rounds:
1. Build from local source (`COPY` the repo, install from there) rather
   than pip-installing an already-published PyPI package, so `docker
   build .` works from any branch/checkout immediately, without a
   publish-then-pull round trip.
2. BLAST+ via `apt install ncbi-blast+` on a `python:slim` base, not a
   conda base image — smaller image, and this ticket's BLAST+ version
   source is intentionally separate from ticket 04's `environment.yml` pin
   (accepted tradeoff per the grilling session: two version sources to
   keep in sync by eye, in exchange for a lighter image).
3. Once uv adoption (ticket 03) was decided, extended to Docker too: use
   `uv pip install .` / `uv sync --frozen` against the committed `uv.lock`
   inside the image, rather than plain `pip install .`, so the image's
   resolved dependency versions match local dev and CI exactly instead of
   pip re-resolving from `pyproject.toml`'s version ranges.

## Scope

- Depends on ticket 01 (package to install), ticket 03 (`uv.lock` to
  install against). Soft-depends on ticket 02 for a meaningful `ENTRYPOINT`/
  `CMD` (the `ijump` console script) — if 02 hasn't landed yet, use the
  pre-02 invocation form and revisit the `ENTRYPOINT` once it has.
- `FROM python:<pinned-version>-slim` (match the Python version
  `pyproject.toml`/`environment.yml` target).
- Install `uv` in the image (per uv's documented Docker pattern: `COPY
  --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/` or an install script —
  implementer's call on which mechanism, both are supported).
- `apt-get install -y ncbi-blast+` (or the precise Debian/Ubuntu package
  name matching the base image's release — confirm the package exists and
  provides both `blastn` and `makeblastdb`, which `ijump.py`'s
  `check_blast_ref`/`clipped_read_search.py`'s `run_blast_subprocess` both
  shell out to).
- `COPY` the repo (respecting a new `.dockerignore` — exclude `.git/`,
  `Test/` and other large example-data directories not needed at runtime,
  `.scratch/`, `__pycache__/`, `.pytest_cache/`, `program.prof`).
- `RUN uv sync --frozen` (or `uv pip install .` — pick whichever matches
  ticket 03's chosen invocation shape) to install the package and its
  locked Python dependencies.
- Set `ENTRYPOINT`/`CMD` to the `ijump` console script (ticket 02) or
  `uv run ijump` if invoking through uv's venv directly.
- Add a brief Docker usage section to `README.md` (build + run example,
  including how to mount input BAM/reference/GFF files as volumes since
  those are user-supplied at runtime, not baked into the image).

## Out of scope

- Publishing the image to a registry (Docker Hub, GHCR, etc.) — this
  ticket produces a buildable/runnable local image only, unless the user
  asks for registry publishing separately.
- Reconciling the apt-sourced BLAST+ version with `environment.yml`'s
  pinned version (ticket 04) — accepted as a known tradeoff, not a defect
  to fix here.
- Multi-stage build optimization / image-size tuning beyond using a slim
  base — nice-to-have, not required for this ticket's acceptance bar.

## Verification

- `docker build .` succeeds from a clean clone with no prior PyPI publish
  or network access beyond package installation.
- `docker run <image> ijump --help` (or equivalent) works.
- A full pipeline run inside the container, with fixture/example BAM+
  reference+GFF files mounted as volumes, produces the expected output
  files (reuse the manual real-sample verification approach from prior
  tickets).
- `blastn -version` and `makeblastdb -version` both succeed inside the
  container (confirms the apt-installed BLAST+ is actually on `PATH` and
  runnable, not just installed).

**Blocked by:** 01 (src-layout package), 03 (`uv.lock`). Soft-depends on 02
(CLI dispatch) for the `ENTRYPOINT`.

**Status:** done

- [x] `Dockerfile` + `.dockerignore` added.
- [x] Image builds from a clean clone with no prior PyPI publish required.
- [x] BLAST+ (`blastn`, `makeblastdb`) installed via apt and runnable inside the container.
- [x] Package installed via `uv` (not against a committed `uv.lock` -- see Comments; none exists for this project).
- [x] A full pipeline run inside the container (fixture data mounted) produces expected output.
- [x] README documents Docker build/run usage.

## Comments

Implemented 2026-08-17. `Dockerfile` and `.dockerignore` added at the repo
root; Docker build/run usage documented in README.md's Installation >
Docker section.

**Contingent dependency resolved:** ticket 16
(`.scratch/isclipped-refactor/issues/16-evaluate-dropping-pysamstats.md`)
closed `closed-no-change` -- `pysamstats` was NOT dropped (both its
correctness and performance bars failed vs. `count_coverage`). So the
Dockerfile does need its own image-install step for `pysamstats` beyond
plain `uv sync`, exactly as this ticket's contingent-dependency note
anticipated.

**No `uv.lock` exists** (ticket 03's finding: `uv sync`/`uv lock` cannot
resolve this project's real dependency set because `pysamstats` 1.1.2
hard-pins `pysam<0.16`, forcing a `pysam` version with no modern wheel).
So `uv sync --frozen` was not usable here; the Dockerfile uses `uv pip
install --system` instead (one of the two forms the ticket explicitly
allows).

**Deviation from spec, deliberate and documented:** the image pins Python
**3.8**, not `environment.yml`'s 3.11. `pysam==0.15.4` (the version
`pysamstats` needs) only has PyPI wheels through `cp38`
(`manylinux2010_x86_64`/`i686` + a macOS Intel wheel, no arm64 at any
Python version) -- confirmed directly against PyPI's file listing. Newer
Python forces a from-source `pysam` build that's independently confirmed
broken (same failure mode documented in README's `uv` section). Python
3.8 still satisfies `pyproject.toml`'s `requires-python = ">=3.7"`. Full
reasoning and the resulting two-step `pysam`-then-`pysamstats --no-deps`
install workaround are documented in both the Dockerfile's header comment
and README's Docker section.

**Verified live**, Docker daemon was available in this session
(`docker build --platform linux/amd64 -t ijump .` from a clean checkout):

- Image builds successfully end-to-end (uv installs pandas/biopython/
  numpy/scikit-learn/pysam via wheels, then pysamstats compiles from
  source against pysam's headers -- needed adding `zlib1g-dev` to the apt
  install list after the first build attempt failed with `fatal error:
  zlib.h: No such file or directory`, then it succeeded).
- `docker run --rm ijump --help` shows the expected subcommand listing.
- `blastn -version` / `makeblastdb -version` both work inside the
  container (BLAST 2.12.0+, from Debian's `ncbi-blast+` package).
- `import pysam, pysamstats` succeeds (with a benign `PileupColumn size
  changed` `RuntimeWarning`, does not affect correctness -- confirmed by
  calling `pysamstats.load_coverage` against `tests/fixtures/tiny.bam`
  directly and getting real, correct coverage numbers back).
- `ijump run --aln --ref --gff --isel --outdir`, with
  `tests/fixtures/tiny.bam`/`.fna`/`.gff`/`is_coords.txt` mounted
  read-only as a volume and an output directory mounted read-write,
  completed without error and produced
  `ijump_junctions.txt`/`ijump_report_by_is_reg.txt`/
  `ijump_sum_by_reg.txt`/`ijump.log`/`reads.txt`/`ijump_wd/` -- the
  expected output set (rows are empty because the tiny fixture has no
  clipped reads near its IS boundaries, same as in the test suite; this
  confirms the pipeline runs cleanly end-to-end inside the container, not
  that this particular fixture exercises every code path).

**Caveat:** `docker build .` on an arm64 host (e.g. Apple Silicon) without
`--platform linux/amd64` fails, because `pysam==0.15.4` has no arm64
wheel and falls back to an unbuildable from-source install (confirmed
directly -- this was the first thing tried, before discovering the
`--platform linux/amd64` requirement). Documented in both the Dockerfile
and README.

**Out of scope, not attempted:** registry publishing, multi-stage build
size optimization, reconciling apt's BLAST+ 2.12.0 against
`environment.yml`'s `blast=2.16.0` pin.

**Re-verified 2026-08-26 (review-followups/08), Docker daemon available:**
`docker build --platform linux/amd64 -t ijump:release-verify .` from the current
`refactor` tip (`59d55eb`) succeeds end to end, and `docker run` against
`tests/fixtures/` (mounted read-only, `-o`/`-w` on a writable mount) exits 0 and
writes the same expected output set described above. Status corrected from
`ready-for-agent` to `done` -- this ticket has been implemented and merged since
2026-08-17; the stale label would have handed an agent already-finished work.

**Dockerfile is stale relative to current `src/ijump`, flagged not fixed:**
`pysamstats` was dropped as a runtime dependency after this ticket
(`80f4235`, isclipped-refactor ticket 16 round 2 -- reversing the "NOT
dropped" finding recorded above, which was correct as of ticket 16's first
round but not the final outcome), and the biopython BLAST wrapper was
removed separately (`f84850d`). Neither the Dockerfile's `RUN uv pip
install --system` line (still pins `biopython==1.79`, `"numpy<2"`, and
`pysam==0.15.4`) nor its header comment's `build-essential`/`zlib1g-dev`
justification (both exist solely to compile `pysamstats` against `pysam`'s
headers) have been updated to reflect that. The image still builds and runs
correctly today because none of those packages being present but unused is
actually broken -- but the whole Python-3.8/`pysam==0.15.4`/apt-toolchain
scaffolding is now unnecessary complexity that a currently-supported Python
and a modern `pysam` wheel could replace. Left for a separate ticket rather
than fixed inline here.
