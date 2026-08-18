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
  Zero hits today: the three `isfinder_parser.py` sites (ticket 01) and the
  `isclipped.py` one (pre-dating this doc, resolved during the isclipped
  refactor) have all been hand-rewritten to `pd.concat`. `README.md` no longer
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
grep -c 'self\.pairs_df' isclipped.py    # 38 — wrong
```

Four of those are `self.pairs_df_path`, a different attribute. ast-grep matching
`self.$ATTR` returns 34, across exactly 3 methods. For a refactor whose whole
question is *which methods share which state*, phantom coupling points you at the
wrong seam. Use ast-grep whenever the answer depends on the code's structure;
grep is still fine for "where is this string".

The same trap in reverse: `$DF.append($$$)` returns 22 hits in this repo, 21 of
them ordinary list appends. Precision comes from writing the *shape* of the idiom
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
' isclipped.py --json=compact
```

Swap `self.$ATTR` for `self.$ATTR = $$$V` to get writers only; the difference
between the two sets is what matters.

### What it found

Of **55** instance attributes, **47 are written once in `__init__` and only read
afterwards** — those are safe to pass as constructor arguments to any extracted
class. Eight are mutated later:

| attribute | readers | written outside `__init__` by |
|---|---|---|
| `min_match` | 4 | `assess_isel_freq`, `report_average` |
| `av_read_len` | 4 | `assess_isel_freq`, `report_average` |
| `blastout_filtered` | 4 | `parseblast` |
| `junctions` | 4 | `call_junctions` |
| `sum_by_region` | 4 | `summary_junctions_by_region` |
| `pairs_df` | 3 | `search_insert_pos` |
| `report_table` | 3 | `report_average` |
| `_index` | 2 | `_crtable_ungapped` |

Six of the eight have a single writer — a clear producer, easy to reason about.
The two with *two* writers are the hazard, and looking at them pays off:

```
isclipped.py:89,91      self.min_match = 150 ; self.av_read_len = 150
isclipped.py:1079-1080  self.min_match = min(self.match_lengths)
                        self.av_read_len = self.read_lengths / self.n_reads_analyzed
isclipped.py:1198-1199  self.min_match = min(self.match_lengths)
                        self.av_read_len = self.read_lengths / self.n_reads_analyzed
```

The computation is duplicated verbatim in `assess_isel_freq` and `report_average`,
each overwriting a `150` placeholder from `__init__`. These are derived values
masquerading as state: whichever method runs first wins, and the placeholder is
live until one of them does. Any extraction that separates those two methods has
to decide who owns the derivation.

## Structural map before reading

The `ast-grep-outline` skill also installed. It is the cheapest way to orient in
a large file:

```bash
ast-grep outline isclipped.py                  # class + its 38 methods
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
