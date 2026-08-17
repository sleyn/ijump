# 04 — Local conda environment + build recipe (not bioconda)

**What to build:** An `environment.yml` for `conda env create` and a local
`meta.yaml` conda-build recipe for `conda build .` / `conda install
--use-local ijump` — reproducible conda installs, not a public bioconda
submission.

**Contingent dependency:** `.scratch/isclipped-refactor/issues/16-evaluate-dropping-pysamstats.md`
is evaluating whether `pysamstats=1.1.2` can be dropped entirely (it has
exactly one call site). Check that ticket's `Status`/`Comments` before
finalizing this ticket's `environment.yml`/`meta.yaml` dependency list —
if 16 lands as a swap, don't pin `pysamstats` here at all.

## Why

Grilled directly with the user: full bioconda submission requires external
review (bioconda-recipes CI/lint + maintainer review, can take weeks, and
an ongoing maintenance commitment) that's out of scope here; a local
recipe gets the reproducibility benefit (BLAST+ pinned alongside Python
deps in one conda-manageable environment) without the publishing overhead.
`README.md`'s existing installation section already leads with a `conda
install` block (lines ~61-75) listing `pysamstats=1.1.2` and friends from
`anaconda`/`bioconda`/`conda-forge` channels — this ticket formalizes that
into reusable files instead of copy-pasted install instructions.

## Scope

- `environment.yml` at the repo root: python version pin, `blast` (BLAST+,
  matching whatever version README's existing conda block specifies —
  re-check `README.md`'s conda section at implementation time), `pysam`,
  `pysamstats=1.1.2` (per README), `pandas`, `biopython`, and any other
  runtime dependency ticket 01's `pyproject.toml` declares. Channels:
  `bioconda`, `conda-forge` (matching README's existing list).
  `conda env create -f environment.yml` should produce a working dev
  environment.
- Local `meta.yaml` (conda-build recipe) — package metadata (name `ijump`,
  version matching `pyproject.toml`), `build` section that installs the
  package (typically `pip install . --no-deps -vv` inside the conda-build
  sandbox, per standard conda-forge/bioconda recipe convention), `requirements`
  section listing the same host/run dependencies as `environment.yml`
  (host: python + build backend from ticket 01; run: blast + the Python
  runtime deps).
- Confirm `conda build .` succeeds locally and `conda install --use-local
  ijump` installs a working `ijump` command (depends on ticket 02's
  console-script existing for this to be meaningful — if 02 hasn't landed,
  verify via the pre-02 invocation form and note the gap).
- Update `README.md`'s installation section to either replace the ad-hoc
  conda install block with `conda env create -f environment.yml`, or note
  both paths (raw install vs. env file) if there's a reason to keep both —
  implementer's call, but don't leave two silently-diverging sets of pinned
  versions (the README block and `environment.yml`) as separate sources of
  truth going forward.

## Out of scope

- Public bioconda submission/PR — explicitly deferred per the grilling
  decision.
- Docker (ticket 05) — the Dockerfile there uses apt for BLAST+, not this
  ticket's conda recipe; they're independent, not layered.
- uv/PyPI publishing (ticket 03) — conda's dependency declarations are
  separate from `uv.lock`; don't try to generate one from the other.

## Verification

- `conda env create -f environment.yml` succeeds from a clean clone.
- `conda build .` succeeds against the local `meta.yaml`.
- `conda install --use-local ijump` (or equivalent) produces a working
  install; run the manual real-sample verification command (reused from
  prior tickets) inside the resulting environment.
- README's install instructions match what's actually in the repo (no
  stale duplicate version pins).

**Blocked by:** 01 (src-layout package, for `pyproject.toml` to exist);
soft-depends on 02 for a meaningful post-install smoke test, same caveat as
ticket 03.

**Status:** done

- [x] `environment.yml` added, pins match README's documented versions (or README updated to match it).
- [x] Local `meta.yaml` conda-build recipe added.
- [x] `conda env create -f environment.yml` succeeds from a clean clone.
- [x] `conda build .` succeeds.
- [x] Manual real-sample run succeeds inside a conda-installed environment.
- [x] README's conda instructions reconciled with `environment.yml` (single source of truth).

## Comments

- Ticket 16 (`pysamstats` drop evaluation) closed `closed-no-change`, so
  `pysamstats=1.1.2` stays pinned in both `environment.yml` and `meta.yaml`
  as instructed.
- `environment.yml`: `python=3.11`, `blast=2.16.0`, `pysam`,
  `pysamstats=1.1.2`, `biopython=1.79` (pinned <1.80 — `Bio.Blast.Applications`
  was removed later), `pandas<3`, `numpy<2`, `scikit-learn`. Channels
  `bioconda`, `conda-forge`, matching README.
- `meta.yaml`: standard noarch:python recipe, `pip install . --no-deps -vv`
  build script, `ijump = ijump.cli:main` entry point (ticket 02), `host`/`run`
  requirements mirroring `environment.yml`, `test: commands: [ijump --help]`.
  Host section pins `python=3.11` explicitly — an unpinned host `python`
  let the solver land on a bleeding-edge Python release with no bioconda
  builds yet for the run-section pins, which produced multi-minute
  backtracking during `conda build .` before eventually still resolving;
  pinning host to match the run environment's Python version avoided that.
- **Real bug found and fixed during manual verification, not just a
  packaging issue:** `region_summary.py`'s
  `temp.at[0, is_coords.keys()] = 0` was never valid pandas usage —
  `.at` only ever supported single-scalar assignment. It happened to
  silently work on ancient pandas (0.x/1.0.x, e.g. the `bioinfo` conda
  env's pandas 1.0.4) but raises `pandas.errors.InvalidIndexError` on
  pandas ≥ roughly 1.1, well before the pandas ≥ 2 cutoff a stale code
  comment in `environment.yml`/`meta.yaml` had previously (incorrectly)
  attributed the failure to. Fixed by switching to `.loc[0,
  list(is_coords.keys())] = 0`, which is valid for multi-column
  assignment on all pandas versions tested (1.5.3 through 3.0.5). This
  let `environment.yml`/`meta.yaml` drop the pandas pin down to `<3`
  instead of staying stuck on `1.5.3` — directly addresses the standing
  "worth switching to a more modern pandas" question. `numpy` is
  separately capped `<2` because `junction_pairing.py` still uses the
  numpy-2.0-removed `np.int0` alias; that's an unrelated, still-open
  modernization item, not touched here (out of scope for a packaging
  ticket — flagging, not fixing, per this ticket-set's convention).
- Manual real-sample run (`ijump run` against `Test/Sample.bam` +
  `Test/A_baumannii_assembly.fna` + the ATCC17978 GFF + ISTable) verified
  successful inside a conda-created `ijump-verify` env (exit 0, full
  `ijump_sum_by_reg.txt`/`ijump_report_by_is_reg.txt` output produced —
  previously it silently truncated to a placeholder-free output at older
  pandas without hitting the `.at` bug, and crashed outright at newer
  pandas before the fix).
- `pytest` in the conda env: 45 passed / 1 failed — the 1 failure
  (`test_read_count_mtx_rejects_invalid_orientation`,
  `ISClipped._read_count_mtx` missing) is the same pre-existing,
  unrelated failure already confirmed by ticket 02's agent.
- `conda build .` succeeded end-to-end (recipe test step ran `ijump
  --help` inside the build sandbox and passed) producing
  `noarch/ijump-1.0.4-py_0.conda`; `conda create -c local -c bioconda -c
  conda-forge ijump` installed it and `ijump --help` worked in the
  resulting environment.
- Environment note: this machine's conda was upgraded mid-ticket from
  4.11.0 (2022, classic solver) to a fresh 26.5.3 install with libmamba
  as the default solver, cutting `conda env create -f environment.yml`
  from 10+ minutes down to ~20-40 seconds. `conda build .` itself still
  took ~20-32 minutes end-to-end even with libmamba — conda-build's
  recipe-render/build-env/host-env/test-env solve steps appear to add
  substantial overhead beyond a plain `conda create` solve on this
  machine; this is a real, reproducible cost of `conda build .`
  specifically, not a config mistake, and worth knowing about before
  relying on it in a tight feedback loop (e.g. CI).
- README's conda section (`## Conda`) documents both `environment.yml`
  (dev environment, then `pip install -e . --no-deps` for the `ijump`
  command) and the `meta.yaml` `conda build .` / `conda install
  --use-local ijump` path as the single reconciled source of truth; no
  stale duplicate ad-hoc `conda install` block remains.
- Follow-up worth a future ticket: `junction_pairing.py`'s `np.int0`
  usage should be updated to `np.intp` (or similar) so `numpy<2` can be
  dropped too, matching the pandas modernization done here.
