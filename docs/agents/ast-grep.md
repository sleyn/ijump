# Structural search with ast-grep

How to use ast-grep in this repo, and what it has already found.

The `ast-grep` skill (installed via `npx skills add ast-grep/agent-skill`) teaches
the general technique — how to write and debug rules. This file carries what is
true *here*: the repo's own rules, the recipes worth rerunning, and the limits.

## The repo's rules

`sgconfig.yml` at the root registers `rules/` and `rule-tests/`. Both are tracked.

```bash
ast-grep scan .     # run the repo's rules
ast-grep test       # check the rules still do what they claim
```

Two rules, each with test cases and a committed snapshot:

- **`pandas-dataframe-append`** — `DataFrame.append` was removed in pandas 2.0.
  Zero hits today: the `isclipped.py` site (pre-dating this doc, resolved during
  the isclipped refactor) was hand-rewritten to `pd.concat`, and the three
  `isfinder_parser.py` sites went away with that module (ticket 02). `README.md` no longer
  pins a pandas version, so a hit here would be a live break, not just an
  upgrade blocker.
- **`identity-compare-literal`** — `is` against a string/number literal. Zero hits
  today; it exists so the `is '+'` strand comparison that lived in a stale `gff.py`
  cannot come back.

Each rule file carries a `note:` explaining why it is written the way it is. Read
that before changing one — both encode a non-obvious detail.

## Why reach for ast-grep over grep

Grep matches text; a Python identifier is not text. Concretely, while mapping
`ISClipped` state:

```bash
grep -c 'self\.pairs_df' src/ijump/isclipped.py
```

Some of those lines are `self.pairs_df_path` — a different attribute that merely
starts with the same characters. ast-grep matching the attribute access itself
does not count them, and reports which methods each hit sits in. For a refactor
whose whole question is *which methods share which state*, phantom coupling
points you at the wrong seam. (Run both; the counts move as the file does, which
is exactly why they are not written down here.)

Use ast-grep whenever the answer depends on the code's structure; grep is still
fine for "where is this string".

The same trap in reverse: `$DF.append($$$)` matches every `.append` in the repo,
nearly all of them ordinary list appends rather than the pandas idiom you were
looking for. Precision comes from writing the *shape* of the idiom
(`$X = $X.append($$$)`), not the name.

## Recipe: the state-coupling matrix

The question the `isclipped.py` refactor keeps needing answered. Metavariables
bound inside an `inside` clause do surface in `--json` output, which is what makes
this work:

```bash
ast-grep scan --inline-rules '
id: touch
language: python
severity: info
rule:
  pattern: self.$ATTR
  inside: {kind: function_definition, stopBy: end, has: {field: name, pattern: $METHOD}}
' src/ijump/isclipped.py --json=compact
```

Swap `self.$ATTR` for `self.$ATTR = $$$V` to get writers only; the difference
between the two sets is what matters.

### What it found

**These counts are a snapshot, and they drift.** The table below was regenerated on
2026-08-25; an earlier version of it named `min_match`, `av_read_len`, `assess_isel_freq`
and `summary_junctions_by_region`, none of which exist any more, and CLAUDE.md points
readers here as though it were current. Re-run the queries above before trusting the
numbers — the technique is the point, not the figures.

Of **30** instance attributes, **13 are written once in `__init__` and only read
afterwards** — those are safe to pass as constructor arguments to any extracted class.
Seventeen are written or mutated later:

| attribute | methods touching it | written outside `__init__` by |
|---|---|---|
| `_depth_by_region` | 2 | `average_depth` |
| `blastout_filtered` | 4 | `run` |
| `cl_read_cov_overlap` | 2 | `run` |
| `clipped_reads` | 3 | `run` |
| `clipped_reads_bwrd` | 3 | `run` |
| `is_clusters` | 5 | `run`, `search_insert_pos` |
| `is_coords` | 4 | `iscollect` |
| `is_table` | 4 | `iscollect` |
| `is_table_fingerprint` | 3 | `run` |
| `junctions` | 4 | `call_junctions` |
| `match_lengths` | 2 | `run` |
| `n_reads_analyzed` | 2 | `run` |
| `pairs_df` | 3 | `run`, `search_insert_pos` |
| `read_lengths` | 2 | `run` |
| `report_table` | 3 | `run` |
| `sum_by_region` | 3 | `run` |
| `unclipped_depth` | 3 | `count_depth_unclipped` |

Fifteen of the seventeen have a single writer — a clear producer, easy to reason about. The
two with *two* writers are the hazard, and looking at them pays off. Both are written by
`run` and again by `search_insert_pos`:

- `is_clusters` is derived from the IS table. `run` builds it up front so a table carrying
  no cluster column stops the run before any work; `search_insert_pos` builds it lazily for
  a caller driving the pairing step on its own. Two writers, and which one fires depends on
  how the object is driven.
- `pairs_df` starts as a placeholder row from `__init__`, is replaced by
  `search_insert_pos`, then rewritten repeatedly by `run` as columns are added. The
  placeholder is live until the first real write — the same shape the retired
  `min_match` case had, where a `150` from `__init__` stood in until whichever method ran
  first overwrote it.

Any extraction separating those methods has to decide who owns the derivation.

## Structural map before reading

The `ast-grep-outline` skill also installed. It is the cheapest way to orient in
a large file:

```bash
ast-grep outline src/ijump/isclipped.py        # the class and its methods
ast-grep outline . --items imports             # dependency direction
```

## Limits

- **Single-file, syntax-only.** No cross-file dataflow, no call graph, no type
  inference. It will not tell you whether `parseblast` writes state that
  `call_junctions` later depends on — that ordering question still needs reading.
- **It finds seams; it does not protect them.** ast-grep's `--update-all` does
  mechanical codemods safely, but nothing it reports tells you whether an
  extraction preserved behaviour. It does not displace the characterization
  harness described in `.scratch/isclipped-refactor/notes/no-test-safety-net.md`,
  which remains ticket #1.
