# iJump image: builds the package from this repo checkout (no dependency on
# a prior PyPI publish), with BLAST+ from apt and the Python dependency
# install driven by uv.
#
# Base Python version: 3.11, matching environment.yml.
#
# `pysamstats` -- the reason a prior version of this image pinned Python 3.8
# and pysam==0.15.4, and built pysamstats from source against apt-installed
# build-essential/zlib1g-dev -- was dropped as a project dependency
# (isclipped-refactor ticket 16 round 2). pysam, pandas, numpy and
# scikit-learn all have prebuilt wheels for this Python/platform combination,
# so `uv sync` resolves and installs them without a C compiler.
FROM python:3.11-slim

# ncbi-blast+ provides both `blastn` and `makeblastdb`, which
# src/ijump/ijump.py's check_blast_ref and
# src/ijump/clipped_read_search.py's run_blast_subprocess shell out to.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ncbi-blast+ \
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

# No uv.lock is committed -- this project doesn't pin one -- so `uv sync
# --frozen` is not usable here; `uv sync` resolves and installs fresh.
RUN uv sync --no-dev --no-editable

ENV PATH="/app/.venv/bin:$PATH"
ENTRYPOINT ["ijump"]
CMD ["--help"]
