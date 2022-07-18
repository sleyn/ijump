# "Average" workflow

In general workflow is following:
1. Find soft-clipped reads near boundaries of mobile elements.
2. Extract unaligned part of reads.
3. BLAST unaligned parts agains reference.
4. Find best hits with >90% identity to the reference with unique highest bitscore. If two ore more hits share the same bit score - the junction considered as ambigious and skipped. Usually this is happening because aligner map reads to several copies of mobile element.
5. Assess frequency of insertion by simple formula:

$\frac{\frac{R_l + R_r}{2 * ( 1 + B_{min} / a_{Rlen})}}{D_t * ( 1 - m_{match} / a_{Rlen} )}$

where

- *$R_l$* - number of reads that support junction to the target on the "left" side of mobile element
- *$R_r$* - number of reads that support junction to the target on the "right" side of mobile element
- *$D_t$* - average depth of caverage of target region
- *$B_{min}$* - minimum length for unaligned parts that were used in the BLAST step
- *$a_{Rlen}$* - average read length
- *$m_{match}$* - minimum length of the read part that could be aligned to reference. Accessed as the minimum of longest clipped part of the read (*e.g.* for read with CIGAR string 10S120M30S *$m_{match}$* is 30).

	$1 + B_{min} / a_{Rlen}$ - correction for reads that were found on the IS element boundary but were not mapped back due to short clipped part.
	
	$1 - m_{match} / a_{Rlen}$ - correction for reads that were not mapped to the IS element boundary but present on the insertion side.

![iJump workflow](https://github.com/sleyn/ijump/img/ijump_workflow.png)

**NOTE:** The "Average" workflow of iJump do not separate different junction events in one target gene. For example, if IS1 element was inserted in gene X three times in different positions, iJump anyway will show one event of unspecified insertion of IS1 element into gene X with summarized frequency.

### Output

#### ijump_junctions.txt

Contains information about junctions for each read. File contains following columns:

* index  
	 order number

* ID  
	 unique identifier

* IS name  
	 mobile element name

* IS pos  
	 what part of the read matches mobile element

* IS chrom  
	 name of contig where mobile element is located in the reference

* Read name  
	 read name where jubction was observed

* Chrom  
	 name of contig where mobile element jumped

* Position  
	 posistion of the junction

* Orientation  
	 orientation of mobile element relative to junction

* Note  
	 mark if junction is in other mobile elements - usually indicates false positive hits

* Locus tag  
	 locus tag of the affected gene; in the case of intergenic region two locus tags will be shown with us_ or ds_ prefixes that indicate upstream or downstream position of the region relative to the genes.

* Gene  
	 trivial name of the affected gene

#### ijump_report_by_is_reg.txt

Long format of frequency estimation. **NOTE: If you have aligner (like BWA-mem) that produses both soft- and hard-clipped reads you should multiply frequency assessments by 2.**

File contains following columns:
* IS Name  
	 mobile element name

* Annotation  
	 locus tag of the affected gene; in the case of intergenic region two locus tags will be shown with us_ or ds_ prefixes that indicate upstream or downstream position of the region relative to the genes.

* Chromosome  
	 name of contig where affected region is located

* Start  
	 start coordinate of affected region

* Stop  
	 end coordinate of affected region

* Frequency  
	 estimated frequency of the mobile element jumps into a genomic region

* Depth  
	 average coverage of the genomic region

If you have several related samples and want to compare them side by side you can copy all *ijump_report_by_is_reg.txt* files in one folder, rename them as *ijump_<*Sample name*>*.txt* and run:

```
python3 combibe_results.py -d [Folder with ijump report files] -o [Output file with the combined table] -g [GFF file. If provided will add functional annotation of the region]
```

This will merge all results in one table.

#### ijump_sum_by_reg.txt

Wide format of frequency estimation. Table shows raw counts of reads that support junctions instead of frequency estimation. **NOTE: If you have aligner (like BWA-mem) that produses both soft- and hard-clipped reads you should multiply frequency assessments by 2.**

* ann  
	 locus tag of the affected gene; in the case of intergenic region two locus tags will be shown with us_ or ds_   

* chrom  
	 name of contig where affected region is located

* start  
	 start coordinate of affected region

* stop  
	 end coordinate of affected region

* mobile element names  
	 raw reads that support junctions

#### CIRCOS files

iJump can create config files (*data* folder) for [CIRCOS](http://circos.ca/) circular diagrams that represent directions of mobile element jumps. Currently commented because of long processing (will be improved soon).

To run CIRCOS you will need to type:
```
circos -config ./data/circos.conf
```

## Simulation test

To assess accuracy of iJump the defined community was simulated. Simulated data mimics several jumps of one of the copy of IS5 element (has several copies in the genome) in Escherichia coli BW25113 genome. The scripts and auxillary files can be found in the **simulation** folder.

The setup of this computational experiment is following:

| Parent | Genome | IS name | IS start | IS stop | Insert position | Tandem | Frequency of insertion | Coverage | Insert repeat |
|--------|--------|---------|----------|---------|-----------------|--------|------------------------|----------|---------------|
| EC_WT  | EC_1   | IS5_10  | 2059640  | 2060834 | 1970720         | 5bp    | 100%                   | 190      | TAAAA         |
| EC_1   | EC_2   | IS5_10  | 2059640  | 2060834 | 1172149         | 5bp    | 25%                    | 250      | GTGCT         |
| EC_1   | EC_3   | IS5_10  | 2059640  | 2060834 | 3846500         | 5bp    | 5%                     | 50       | CACCG         |
| EC_1   | EC_4   | IS5_10  | 2059640  | 2060834 | 240955          | 5bp    | 1%                     | 10       | GTCGC         |
| EC_1   | EC_5   | IS5_10  | 2059640  | 2060834 | 3358463         | 5bp    | 50%                    | 500      | GCAAT         |

Reads were simulated with [ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)

Alignment was made with [bwa-mem](http://bio-bwa.sourceforge.net/)

BAM file manipulations were performed with [samtools](http://samtools.sourceforge.net/)

Results are following:

| IS Name | Annotation   | Chromosome | Start   | Stop    | Frequency | Depth       |
|---------|--------------|------------|---------|---------|-----------|-------------|
| IS5_9   | BW25113_1117 | CP009273.1 | 1172083 | 1172777 | 2.94%     | 976         |
| IS5_10  | BW25113_1117 | CP009273.1 | 1172083 | 1172777 | 20.90%    | 976         |
| IS5_9   | BW25113_3216 | CP009273.1 | 3356166 | 3358544 | 6.17%     | 940         |
| IS5_7   | BW25113_3216 | CP009273.1 | 3356166 | 3358544 | 0.07%     | 940         |
| IS5_10  | BW25113_3216 | CP009273.1 | 3356166 | 3358544 | 41.39%    | 940         |
| IS5_9   | BW25113_0223 | CP009273.1 | 240814  | 241552  | 0.14%     | 966         |
| IS5_10  | BW25113_0223 | CP009273.1 | 240814  | 241552  | 0.48%     | 966         |
| IS5_1   | BW25113_1890 | CP009273.1 | 1970513 | 1971397 | 0.07%     | 900         |
| IS5_8   | BW25113_1890 | CP009273.1 | 1970513 | 1971397 | 0.07%     | 900         |
| IS5_5   | BW25113_1890 | CP009273.1 | 1970513 | 1971397 | 0.07%     | 900         |
| IS5_4   | BW25113_1890 | CP009273.1 | 1970513 | 1971397 | 12.29%    | 900         |
| IS5_10  | BW25113_1890 | CP009273.1 | 1970513 | 1971397 | 83.29%    | 900         |
| IS5_10  | BW25113_3671 | CP009273.1 | 3844456 | 3846145 | 4.35%     | 966         |

We see that some reads were aligned to other copies of IS5 element. However majority of reads were aligned correctly.
If we summarize all frequences for each affected gene we will get results close to the expected:

| Gene         | Start   | Stop    | Observed | Expected |
|--------------|---------|---------|----------|----------|
| BW25113_1890 | 1970513 | 1971397 | 95.79%   | 100.00%  |
| BW25113_3216 | 3356166 | 3358544 | 47.63%   | 50.00%   |
| BW25113_1117 | 1172083 | 1172777 | 23.84%   | 25.00%   |
| BW25113_3671 | 3844456 | 3846145 | 4.35%    | 5.00%    |
| BW25113_0223 | 240814  | 241552  | 0.62%    | 1.00%    |

Simulation shows that iJump tend to slightly decrease frequency. This is happening because instead of work with individual junction sites it summarizes all sites along the genetic element where junctions were found and assess frequency with average depth in the element. However near junctions coverage has slight drop. The reason of this effect is that aligners have a limit on a length of aligned part of a read. iJump estimates this limit and introduces correction cefficients.
