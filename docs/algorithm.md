# iJump — Algorithm Description

## Table of Contents

- [Overview](#overview)
- [Inputs](#inputs)
- [Shared Steps](#shared-steps)
  - [Step 1 — IS Boundary Windows](#step-1--is-boundary-windows)
  - [Step 2 — Forward Clipped Read Collection (IS→Ref)](#step-2--forward-clipped-read-collection-isref)
  - [Step 3 — BLAST Search (IS→Ref)](#step-3--blast-search-isref)
  - [Step 4 — Junction Table](#step-4--junction-table)
- [Average Mode](#average-mode)
  - [Step 5A — Summarize Junctions by Genomic Region](#step-5a--summarize-junctions-by-genomic-region)
  - [Step 6A — Frequency Estimation](#step-6a--frequency-estimation)
- [Precise Mode](#precise-mode)
  - [Step 5B — Cluster Reference Positions](#step-5b--cluster-reference-positions)
  - [Step 6B — Backward Clipped Read Collection (Ref→IS)](#step-6b--backward-clipped-read-collection-refis)
  - [Step 7B — BLAST Search (Ref→IS)](#step-7b--blast-search-refis)
  - [Step 8B — Junction Pairing](#step-8b--junction-pairing)
  - [Step 9B — Depth Counting at Junction Positions](#step-9b--depth-counting-at-junction-positions)
  - [Step 10B — Frequency Estimation](#step-10b--frequency-estimation)
- [Multi-sample Comparison](#multi-sample-comparison-combine_resultspy)
- [Python Library Usage](#python-library-usage)
  - [scikit-learn](#scikit-learn)
  - [NumPy](#numpy)
  - [pandas](#pandas)

---

## Overview

iJump detects Insertion Sequence (IS) element rearrangements in evolved bacterial populations from short-read (Illumina) sequencing data. The core signal it exploits is **soft-clipped reads**: when an IS element has jumped into a new genomic location, reads that span the insertion junction will be partially aligned — one part maps to the reference at the insertion site, and the unaligned (clipped) part maps to the IS element itself.

iJump provides two workflows:

- **Average mode** — estimates frequency at the gene/intergenic-region level; more sensitive, less precise.
- **Precise mode** — localises each insertion to exact base-pair coordinates before frequency estimation; more accurate, more computationally intensive.

Both workflows share the same initial steps for clipped read collection and BLAST alignment.

---

## Inputs

| Input | Format | Description |
|---|---|---|
| BAM file | BAM/SAM | Short reads aligned to the reference genome |
| Reference genome | FASTA | One record per contig |
| Genome annotations | GFF3 | PATRIC/PROKKA-style, with `##sequence-region` headers |
| IS element coordinates | TSV | `IS_name  contig  start  stop` |

---

## Shared Steps

### Step 1 — IS Boundary Windows

*Implemented in `isclipped.py` → `set_is_boundaries()`*

For each IS element, a search window of `radius` bp (default 200 bp) is placed symmetrically around both its **start** and **stop** boundaries in the reference. If a window would extend past the end of a contig, it is split into two intervals that wrap around, treating the contig as circular.

```
[start − radius,  start + radius]   ← start-side window
[stop  − radius,  stop  + radius]   ← stop-side window
```

### Step 2 — Forward Clipped Read Collection (IS→Ref)

*Implemented in `isclipped.py` → `collect_clipped_reads()` → `_crtable_ungapped()`*

All reads overlapping each search window are fetched from the BAM file. A read is kept only if:

1. It is mapped.
2. Its CIGAR string contains `S` (soft-clip).
3. The clipped segment is on the **correct side** of the IS element:
   - Reads near the **start** boundary must have a **left-clipped** segment (unaligned part precedes the aligned part in the read).
   - Reads near the **stop** boundary must have a **right-clipped** segment (unaligned part follows the aligned part).

For each qualifying clipped segment the following is recorded: clipped sequence, reference coordinate of the junction nucleotide, side of clipping, read orientation, and the associated IS element. Read lengths and the lengths of matched (M) CIGAR segments are accumulated for later correction coefficients.

### Step 3 — BLAST Search (IS→Ref)

*Implemented in `isclipped.py` → `runblast()` and `parseblast()`*

All clipped sequences of length ≥ 10 bp are written to FASTA and searched against the reference genome with **BLASTn** (e-value ≤ 0.001, word size 10, tabular output format 6).

Hits are filtered:

| Filter | Criterion |
|---|---|
| Minimum identity | ≥ 75% |
| Best hit selection | Only the hit(s) with maximum bitscore per query read are kept |
| Uniqueness (IS→Ref only) | If two or more hits share the top bitscore, the read is discarded as ambiguous (typically caused by IS element copies) |

The reference position of the junction (left or right end of the BLAST hit, depending on clip orientation) is stored for each passing read.

### Step 4 — Junction Table

*Implemented in `isclipped.py` → `call_junctions()`*

Each filtered BLAST hit becomes one row in the **junction table** with fields:

- IS element name and its contig
- Reference contig and position of the junction
- Orientation (whether the non-clipped part of the read is to the left or right of the junction)
- GFF annotation at that position (locus tag, gene name)
- A flag marking junctions that map back into IS element boundaries (these are within-IS artefacts and are excluded from frequency estimation)

---

## Average Mode

### Step 5A — Summarize Junctions by Genomic Region

*Implemented in `isclipped.py` → `summary_junctions_by_region()`*

Junctions flagged as within-IS artefacts are removed. For each remaining junction, the GFF annotation it falls within is looked up. Supporting read counts are accumulated in a wide-format summary table indexed by genomic region, with one column per IS element.

### Step 6A — Frequency Estimation

*Implemented in `isclipped.py` → `report_average()`*

For each (IS element, genomic region) pair with at least one supporting read:

$$
\text{Frequency} = \frac{\dfrac{R_l + R_r}{2} \cdot \left(1 + \dfrac{B_{\min}}{a_{rlen}}\right)}{D_t \cdot \left(1 - \dfrac{m_{match}}{a_{rlen}}\right)}
$$

| Symbol | Meaning |
|---|---|
| $R_l + R_r$ | Total clipped reads supporting junctions in the region (left + right edges of IS) |
| Division by 2 | Each insertion creates two junctions; dividing by 2 avoids double-counting |
| $B_{\min}$ | Minimum clipped length accepted for BLAST (10 bp); corrects for excluded very-short clips |
| $a_{rlen}$ | Average read length |
| $D_t$ | Average coverage depth of the genomic region (computed with pysamstats) |
| $m_{match}$ | Minimum observed matched-segment length; corrects for reads too short to produce a detectable clip |

The correction factors account for two systematic biases:

- $\left(1 + B_{\min}/a_{rlen}\right)$ — adds back reads whose clipped part was shorter than the BLAST minimum and were therefore excluded.
- $\left(1 - m_{match}/a_{rlen}\right)$ — removes from the denominator the fraction of reads at the insertion site that cannot produce a detectable clipped segment because their matched part is too short for the aligner.

**Output**: `ijump_report_by_is_reg.txt`

> **Note:** Average mode does not separate multiple insertion events within the same gene. If IS1 inserted into gene X at three different positions, iJump reports one combined event for IS1 in gene X with the summed frequency.

---

## Precise Mode

### Step 5B — Cluster Reference Positions

*Implemented in `isclipped.py` → `make_gene_side_regions()`*

BLAST hits from Step 3 that fall within IS element boundaries are removed. The remaining hit positions represent putative insertion sites in the genome. These positions are clustered independently per chromosome using **hierarchical agglomerative clustering** (single linkage, distance threshold = 30 bp). Each resulting cluster defines a compact **reference region** (min − 5 bp to max + 5 bp) likely to contain one or more IS insertions.

### Step 6B — Backward Clipped Read Collection (Ref→IS)

*Implemented in `isclipped.py` → `crtable_bwds()`*

Within each reference region from Step 5B, all soft-clipped reads are collected again (without the IS-side directionality filter). The clipped parts of these reads should map to IS elements. In addition, the coverage contributed by the **aligned** portions of these clipped reads is tallied per reference position — this "overlap coverage" is used later to avoid double-counting.

### Step 7B — BLAST Search (Ref→IS)

*Implemented in `isclipped.py` → `runblast()` and `parseblast()`*

The same BLAST pipeline as Step 3 is applied to the backward reads. Here the uniqueness filter is relaxed: when multiple hits share the top bitscore, the first one is kept rather than discarding the read, because the IS element assignment is handled proportionally in Step 10B.

### Step 8B — Junction Pairing

*Implemented in `isclipped.py` → `search_insert_pos()` → `_find_pair()`*

IS element copy suffixes (e.g. `IS1_1`, `IS1_2`) are stripped so all copies of the same IS family are grouped together. For each (IS element family, chromosome) group, left and right junction positions are paired to identify the two edges of a single insertion event:

1. A **closeness matrix** $C$ is constructed where $C_{ij} = 1$ if left position $i$ and right position $j$ are within `max_is_dup_len` = 20 bp of each other (the expected target-site duplication length). Wrap-around proximity near contig ends is also checked.
2. Positions are grouped into clusters based on overlapping proximity columns.
3. Within each cluster, positions are sorted by supporting read count (descending).
4. Positions are greedily paired: each left junction is matched to the closest-count right junction from the same cluster (penalising cross-cluster matches by 10 000).
5. Left or right junctions with no partner become **orphan** observations with the missing coordinate set to 0.

### Step 9B — Depth Counting at Junction Positions

*Implemented in `isclipped.py` → `count_depth_unclipped()`*

For each junction position, reads with no soft- or hard-clipping in their CIGAR string are counted. Reads where the junction falls too close to the read end (within `Dist`, the distance between left and right junctions) are excluded because it is impossible to determine whether they are truly non-clipped or merely have their clip at the read edge.

### Step 10B — Frequency Estimation

*Implemented in `isclipped.py` → `assess_isel_freq()`*

**Read attribution** — When multiple IS elements share a junction position, the raw read count at that position is split proportionally among IS elements according to their backward-BLAST read counts:

$$
N_{IS_i, pos} = N_{raw, pos} \cdot \frac{C_{IS_i, pos}}{\sum_j C_{IS_j, pos}}
$$

**Read count correction** — All clipped read counts are inflated to compensate for reads whose clipped part is too short to be detected:

$$
N_{cl,corrected} = \frac{N_{cl}}{1 - m_{match}/a_{rlen}}
$$

**Overlap subtraction** — At the left junction, reads clipped at the right junction that extend across the left junction position are subtracted from the overlap coverage (and vice versa), to avoid counting the same read twice when the two junctions are close:

$$
N_{ov,formula,l} = \begin{cases}
N_{ov,l,corrected} - N_{cl,r,corrected} & \text{if } pos_r > pos_l > 0 \\
N_{ov,l,corrected} & \text{otherwise}
\end{cases}
$$

**Per-junction frequency**:

$$
Freq_l = \frac{N_{cl,l,corrected}}{N_{ncl,l} + N_{ov,formula,l} + N_{cl,l,corrected} + 0.1}
$$

| Symbol | Meaning |
|---|---|
| $N_{cl,corrected}$ | Corrected clipped read count supporting this junction for the given IS element |
| $N_{ov,formula}$ | Overlap clipped reads from the opposite junction, after correction and subtraction |
| $N_{ncl}$ | Unclipped reads at the junction position (background depth) |
| 0.1 | Pseudocount to avoid division by zero |

**Final frequency** is the average of left and right junction frequencies when both are present, or whichever is available for orphan observations:

$$
Freq = \begin{cases}
Freq_r & \text{if } Freq_l = 0 \\
Freq_l & \text{if } Freq_r = 0 \\
(Freq_l + Freq_r) / 2 & \text{otherwise}
\end{cases}
$$

**Output**: `ijump_junction_pairs.txt`

---

## Multi-sample Comparison (`combine_results.py`)

When multiple samples from the same experiment are available, `combine_results.py` merges their individual iJump outputs into a single comparative table:

1. **Collect** all `ijump_*.txt` report files from a directory.
2. **Filter** low-confidence observations (depth ≤ 10 reads).
3. **Outer-join** all per-sample DataFrames on the identity columns (IS name, chromosome, start, stop) so every observed IS-region pair appears as one row, with per-sample frequencies as columns and 0 for absent observations.
4. **Dense mode** (precise only): orphan single-edge observations are merged with complete pairs that share the same known junction coordinate, then read counts are summed.
5. **Annotate** with GFF (locus tags, gene names, functional products), looking up the annotation at the midpoint (average mode) or at each junction independently (precise mode).
6. **Classify** each insertion as:
   - **A** (acquired) — not observed in unevolved ancestor samples
   - **P** (pre-existing) — present in unevolved ancestor samples
7. **Filter** to events with MAX frequency ≥ 1 % for the "selected" output table.

---

## Python Library Usage

### scikit-learn

Only one function from scikit-learn is used, applied in a single targeted step:

**`sklearn.cluster.AgglomerativeClustering`** — `isclipped.py`, `_hclust()` / `make_gene_side_regions()`

```python
AgglomerativeClustering(n_clusters=None, distance_threshold=30, linkage='single')
    .fit(X.to_numpy().reshape(-1, 1))
```

- **Input**: 1D array of genomic positions on one chromosome, reshaped to a column vector for sklearn.
- **Purpose**: groups junction candidates within 30 bp of each other (single-linkage distance) to define compact insertion-site regions for the backward search pass (Step 5B).
- **Called via**: `groupby('Chrom')['Position'].transform(_hclust)` — applied independently per chromosome inside a pandas groupby.

---

### NumPy

NumPy is used in `isclipped.py` for all array-level numerical operations, primarily in the junction-pairing and frequency-estimation logic.

| Operation | Location | Purpose |
|---|---|---|
| `np.zeros(shape)` | `_find_pair`, `_read_count_mtx` | Initialise closeness matrix and per-IS read-count matrices |
| `np.ones_like(pos_r)` | `_find_pair:715` | Broadcast a comparison against all right positions simultaneously |
| `np.abs()` | `_find_pair`, `count_depth_unclipped` | Absolute distance between positions |
| `np.argmin()` | `_find_pair:779` | Greedy selection of right junction with the closest read count to the left junction |
| `np.argsort()` | `_find_pair:756` | Sort positions within a cluster by read count (descending) |
| `np.where()` | `_find_pair` | Index positions and counts belonging to a given cluster |
| `np.unique()` | `_find_pair:752` | Enumerate distinct cluster IDs |
| `np.any()` / `np.all()` | `keep_pair`, `_find_pair:738` | Vectorised membership tests: check if a position falls in any region interval; check if a closeness matrix column is non-zero |
| `np.sum()` | `_find_pair:721,776` | Count occupied cells in closeness matrix |
| `np.arange()` | `_find_pair:772` | Create index array to track orphan right-junction positions |
| `np.array()` | `fisher_test_clr_number`, `count_depth_unclipped` | Build Fisher contingency table; collect read-edge positions |
| `np.min()` | `count_depth_unclipped:956` | Distance from junction to nearest read edge |
| `.reshape(-1, 1)` | `_hclust`, `_resore_orig_counts` | Format 1D array for sklearn input; row-normalise 2D count matrix |
| Row normalisation `mtx /= mtx.sum(1).reshape(-1,1)` | `_resore_orig_counts:927` | Convert raw IS-attributed counts to proportions for read-attribution |
| `np.int0` cast | `_find_pair:715` | Ensure closeness matrix stores integer flags |

---

### pandas

pandas is the universal data layer. Every internal table (clipped reads, junctions, pairs, summary, report) is a DataFrame, and all inter-step data transfer uses DataFrame operations.

#### Table construction

| Operation | Location | Purpose |
|---|---|---|
| `pd.DataFrame(columns=[...])` | `isclipped.py` init methods | Declare typed empty tables for clipped reads, junctions, pairs, report |
| `pd.DataFrame.from_dict(..., "index")` | `ijump.py:187` | Build clipped-reads table from a dictionary of rows — faster than row-by-row append for unknown row counts |
| `pd.concat([...])` | `ijump.py:264`, `search_insert_pos` | Merge left and right position frames for depth counting; concatenate per-IS pair tables into final pairs frame |

#### Filtering and selection

| Operation | Location | Purpose |
|---|---|---|
| `.query(expr)` | Throughout | Filter by orientation (`"left"/"right"`), note (`"!= 'IS element'"`), position (`"> 0"`), depth (`"> 10"`), frequency (`">= 0.01"`) |
| `groupby(...)[col].transform(max) == col` | `parseblast:508` | Vectorised best-hit selection: mark rows whose bitscore equals the per-query maximum |
| `.drop_duplicates()` | `ijump.py`, `isclipped.py` | Remove redundant position entries before depth counting |
| `.drop(columns=[...])` / `.rename(columns={...})` | Throughout | Column normalisation before joins |

#### Aggregation and reshaping

| Operation | Location | Purpose |
|---|---|---|
| `.groupby().aggregate(['min','max'])` | `make_gene_side_regions:588` | Extract region boundaries (min/max position) from clustered junction points |
| `.groupby().transform(_hclust)` | `make_gene_side_regions:583` | Apply hierarchical clustering independently per chromosome |
| `.groupby().count().reset_index()` | `search_insert_pos:829` | Count reads per unique (position, IS, orientation) combination |
| `pd.melt(id_vars=..., var_name=..., value_name=...)` | `report_average:1200` | Pivot the wide sum-by-region table (IS elements as columns) to long format for frequency calculation |
| `.groupby().agg('sum')` | `combine_results.py:349` | Collapse IS element copies (IS1_1 + IS1_2 → IS1) by summing frequencies |
| `reduce(lambda df1, df2: pd.merge(..., how='outer'), dfs)` | `combine_results.py:107` | Iterative outer-join of all per-sample report DataFrames into one wide comparative table |

#### Row-wise computation

| Operation | Location | Purpose |
|---|---|---|
| `.apply(lambda row: ..., axis=1)` | `assess_isel_freq`, `report_average`, `filter_pairs` | Frequency formula evaluation, annotation lookup, region membership test |
| `.itertuples(index=False)` | `call_junctions`, `assess_isel_freq` | Efficient row iteration when building junction table and reading count matrices |
| `.apply(max, axis=1)` | `combine_results.py:261` | Per-row maximum frequency across samples |

#### I/O

| Operation | Location | Purpose |
|---|---|---|
| `pd.read_csv(path, sep='\t')` | `parseblast`, `combine_results.py` | Read BLAST tabular output and per-sample iJump reports |
| `.to_csv(path, sep='\t', index=False)` | `ijump.py`, `isclipped.py` | Write all result tables (junctions, pairs, sum-by-region, report) as TSV |
| `.fillna(0)` | `combine_results.py:110` | Fill zeros for samples that did not observe a given insertion |
