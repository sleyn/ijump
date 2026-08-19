# 02 — Fix frequency_estimation.py's positional row-indexing under pandas 3.0

Status: done

**What to build:** `src/ijump/frequency_estimation.py:234-235` does

```python
pairs_df["Frequency"] = pairs_df[["Frequency_l", "Frequency_r"]].apply(
    lambda x: _calc_freq_precise(x[0], x[1]), axis=1
)
```

Inside the lambda, `x` is a row `Series` indexed by column *labels*
(`"Frequency_l"`, `"Frequency_r"`), not position. Older pandas fell back to
positional indexing for `x[0]`/`x[1]` (with a `FutureWarning`); pandas 3.0
removed that fallback, so `Series.__getitem__` does a pure label lookup and
raises `KeyError: 0` (there is no column literally named `0`). `pyproject.toml`
pins plain `"pandas"` with no ceiling, so a fresh install already resolves to
pandas 3.0.5 and hits this — `estimate_frequencies` (PRECISE mode's
frequency calc) is broken there today, not just when someone upgrades later.

Fix `frequency_estimation.py` to index by label (`x["Frequency_l"],
x["Frequency_r"]`) or by explicit position (`x.iloc[0], x.iloc[1]`) instead
of bare integer `__getitem__`, so the call is correct on every currently
supported pandas version, not just the older ones with the positional
fallback.

**Blocked by:** None — can start immediately

- [x] `tests/test_estimate_frequencies.py::test_estimate_frequencies_matches_pinned_golden_output`
      and `::test_estimate_frequencies_does_not_mutate_input_pairs_df` pass
      against a freshly resolved (unpinned) pandas install, not just an
      environment with an old pandas pinned by hand
- [x] No behavior change versus the pinned golden output on a pandas version
      where the old code did work (e.g. the conda-resolved pandas in
      `environment.yml`, if older than 3.0) — this is a compatibility fix,
      not a logic change
- [x] `ruff check`/`mypy` on `src/ijump` and `tests` still pass

## Comments

Switched the lambda to label-based indexing (`x["Frequency_l"], x["Frequency_r"]`)
instead of bare integer `__getitem__` — correct under both the old positional-fallback
behavior and pandas 3's label-only `Series.__getitem__`, so no version-conditional code
needed. Verified `tests/test_estimate_frequencies.py` (both cases) against a fresh
pip-installed pandas 2.3.3 (the newest currently on PyPI — pandas 3.0 isn't published
yet, so the ticket's "already resolves to pandas 3.0.5" couldn't be reproduced verbatim
today, but label lookup is unaffected by which side of that release the fix runs on).
Since the indexing bug is now fixed for any pandas version, also dropped the now-stale
`pandas<3` pin and its comment from `environment.yml`.
