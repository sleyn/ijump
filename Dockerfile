# iJump image: builds the package from this repo checkout (no dependency on
# a prior PyPI publish), with BLAST+ from apt and the Python dependency
# install driven by uv.
#
# Base Python version: 3.8, NOT environment.yml's 3.11.
#
# Why: pysamstats 1.1.2 (PyPI, last released 2018, unmaintained -- see the
# "uv" subsection of README.md's Installation > Docker section for the full
# writeup) hard-pins `pysam<0.16` in its own package metadata. Modern pysam
# has no wheel compatible with that old pin, and the pysam 0.15.x sdist
# fails to build on current toolchains/Python (verified independently in
# ticket 03: fails on `pkg_resources`, then a deeper Cython/sdist gap, even
# with build tooling pre-installed). pysam 0.15.4's *last* PyPI wheel
# lineage tops out at cp38 (manylinux2010_x86_64) -- confirmed directly
# against PyPI's file listing for that release. Pinning this image to
# Python 3.8 (still satisfies pyproject.toml's `requires-python = ">=3.7"`)
# is what makes `pysam==0.15.4` installable from a prebuilt wheel instead of
# a from-source build, which in turn is what makes pysamstats' own
# source-only sdist (it ships a precompiled opt.c, so no Cython needed, but
# it does need pysam's C headers to match at compile time) buildable with
# just a C compiler. This is the "own image-install step beyond plain uv
# sync" the ticket anticipated needing.
#
# Verified end-to-end with a live `docker build --platform linux/amd64 .`
# and `docker run`: image builds, `blastn -version`/`makeblastdb -version`
# both work, `import pysam, pysamstats` succeeds (with a benign
# `PileupColumn size changed` RuntimeWarning -- pysamstats' opt.c was
# compiled against pysam 0.15.4's struct layout, which differs slightly
# from what its .pyx expects at the Python-object level; it does not
# affect correctness, confirmed by calling pysamstats.load_coverage
# against tests/fixtures/tiny.bam and getting real, correct-looking
# coverage numbers back), and `ijump run` completes against
# tests/fixtures/tiny.bam/.fna/.gff mounted as volumes, producing the
# expected output files. IMPORTANT CAVEAT: this only works on
# linux/amd64. pysam 0.15.4's PyPI wheels were built for a 32/64-bit
# x86 world (cp38-manylinux2010_{x86_64,i686} and a macOS Intel wheel)
# and there is no arm64 wheel at any Python version for that release, so
# `docker build .` on an Apple Silicon host (or any arm64 builder)
# without `--platform linux/amd64` falls back to a from-source pysam
# build, which fails (no cython, no precompiled .c in the sdist) exactly
# like the failure already documented in README's uv section. Pass
# `--platform linux/amd64` explicitly when building on arm64.
FROM python:3.8-slim

# ncbi-blast+ provides both `blastn` and `makeblastdb`, which
# src/ijump/ijump.py's check_blast_ref and
# src/ijump/clipped_read_search.py's run_blast_subprocess shell out to.
# build-essential is needed to compile pysamstats' bundled opt.c against
# pysam's headers (see above); it is not removed afterwards -- multi-stage
# / image-size tuning is explicitly out of scope for this ticket.
# zlib1g-dev supplies zlib.h, which pysam's bundled htslib headers
# (transitively included by pysamstats/opt.c) need at compile time --
# without it the pysamstats build fails with `fatal error: zlib.h: No
# such file or directory` (hit and fixed during verification).
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ncbi-blast+ \
        build-essential \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install uv itself (astral's documented Docker pattern) rather than
# `pip install uv`, so the image always gets the tool's own maintained
# static build.
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app

# Copy the repo checkout in (respecting .dockerignore) and install from it,
# so `docker build .` works from any branch/checkout without a
# publish-then-pull round trip through PyPI.
COPY . .

# No uv.lock is committed (see README's uv section: `uv sync`/`uv lock`
# cannot resolve this project's real dependency set today, precisely
# because of the pysam/pysamstats situation above), so `uv sync --frozen`
# is not usable here. Install the resolvable dependencies with `uv pip
# install --system` instead, then install pysam/pysamstats as a deliberate
# two-step workaround for the pin conflict described above, then install
# the `ijump` project itself with --no-deps since every dependency it
# declares has already been installed explicitly above.
RUN uv pip install --system \
        pandas \
        "biopython==1.79" \
        "numpy<2" \
        scikit-learn \
        "pysam==0.15.4" \
    && uv pip install --system --no-deps "pysamstats==1.1.2" \
    && uv pip install --system --no-deps .

ENTRYPOINT ["ijump"]
CMD ["--help"]
