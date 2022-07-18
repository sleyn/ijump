# iJump

<img src="./img/logo/logo1.png" width="150" height="150">

Software for search of Insertion Sequences (IS) rearrangements in evolved population sequencing data.
*Program is in a development stage, and I will appreciate any feedback.*

**iJump** searches for IS elements rearrangements in evolved populations of single organism and estimates what fraction of a population is affected by the rearrangement. iJump uses soft-clipped reads to find evidense for rearrangements. 

**NOTE:** Working with short-read-only assembled genomes is difficult with iJump. The reason is that usually IS elements are repetitive regions which are difficult to resolve for assemblers. This often result in shreading IS elements to several/many sometimes overlapped short contigs. This introduces difficulty either for boundaries determination and for mapping algorithms.

## Motivation

Genome rearrangements and especially IS rearrangements are powerful tools for evolution in all domains of life. Currently very few software tools exists that can work with data of mixed populations, not the clonal ones.
The goal of the iJump is to estimate frequency of rearrangements in the evolved population. The scenario that was in our mind is following:

1) Bacterial population is going through experimental evolution.
2) During the evolution process the initial population split into several subpopulations.
3) In subpopulations various IS rearrangements occur.

The only tool for bacterial population data that was found is a [*breseq*](http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing). However although it is an excellent and comprehensive tool we've found couple limitations of its use. (1) It very slow on high coverage data and (2) it is impossible to use only selected analysis tools.

## Installation

**iJump** do not need special installation. Just clone repository:
```
git clone https://github.com/sleyn/ijump.git
```

But it is dependent on several Python libraries:
* **biopython**
* **pandas**
* **pysam**
* **pysamstats**
* **numpy**
* **scipy**
* **sklearn**

### Conda

The dependences could be installed in Conda environment:
```
conda install \
    -c anaconda \
    -c bioconda \
    -c conda-forge \
    biopython=1.79 \
    pandas=1.3.5 \
    numpy=1.21.0 \
    scikit-learn=0.24.2 \
    pysam=0.15.3 \
    pysamstats=1.1.2 \
    scipy=1.4.1
```

### Docker

**----IN DEVELOPMENT----**

 ## Usage

 ### Input

 iJump requires four files for input:
 1. File with mobile elements coordinates
 2. Reference DNA contigs fasta file.
 3. GFF file with reference genome annotations.
 4. BAM file of aligned Illumina reads.

 ![iJump input and output](https://github.com/sleyn/ijump/img/ijump_input.png)

 #### Mobile elements coordinates file

 File with mobile elements coordinates shoud be tab-separated tables of the following structure:
 ```
 IS_Name    Contig_Name Start_position  End_position
 ```

 For example:
 ```
 ISAcsp3	NODE_1	2980551	2981283
 ```

If you don't have file with coordinates of mobile elements you can:

1) Preferred. Do manual BLAST against standalone ISFinder database. Database could be downloaded from:

- [Author GitHub](https://github.com/thanhleviet/ISfinder-sequences)
- [My Fork](https://github.com/sleyn/ISfinder-sequences) with already built BLASTn database.

Do BLASTn search:

```
blastn -query <Genome> -db <BLASTn database from IS.fna> -out <Output file> -outfmt 6
```

Parse the output table with **isfinder_db_parcer.py** script:
```
python3 isfinder_db_parcer.py -b <BLAST output in outfmt 6 format> -o <Output directory>
```


2) Find them from ISFinder website using their [BLAST](https://isfinder.biotoul.fr/blast.php) against your reference contigs.

It will return you  html page of hits that you can download and parse with **isfinder_parser.py**:
```
python3 isfinder_parse.py -i <ISfinder BLAST HTML page>
```

Both parsers will find non-overlapping hits with empirical E-value threshold 1E-30.

**NOTE:** It was observed that if the contig FASTA header (the line starting with ">") is long then ISFinder BLAST does not produce "Query=" string with the contig name.  This line is critical for `isfinder_parse.py` work. If the script reports empty table please change header by using sorter contig names or removing auxiliary information.

#### Reference Fasta

Regular Fasta file with one Fasta-record per contig:

```
>Contig1
gctagctagctagctacgtagctagctagctacgtacgtacgtagcta...

>Contig2
cgtagctgctagctagctagctagcgtacgtacgtagctacgtacgta...

...
```

#### GFF file

iJump is working with it's own gff module that is tuned for PATRIC/PROKKA-style GFF.

Some specifications:
- The comment string `##sequence-region	accn|[contig name]	[contig start coordinate]	[contig end coordinate]` is required for each contig in the reference.
- Info field should have the following fields:
  	
	- ID
	- product
	- locus_tag
	- (optional, if gene has trivial name) gene

Example:

```
##gff-version 3								
##sequence-region	accn|NODE_1_length_3909467_cov_533.478_ID_22129	1	3909467					
NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	34	336	.	-	1	ID=fig|400667.82.peg.1;product=hypothetical protein;locus_tag=AUO97b_00141
NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	352	1578	.	-	1	ID=fig|400667.82.peg.2;product=phage replication protein Cri;locus_tag=AUO97b_00142;gene=cri
NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	1724	2098	.	+	2	ID=fig|400667.82.peg.3;product=helix-turn-helix family protein;locus_tag=AUO97b_00143
```

If you have another style of GFF unfortunately you have to reformat it at the current stage of development.

#### BAM file

BAM file with aligned short reads. BAM file should contain soft-slipped reads. On the current stage iJump don't use hard-clipped reads.

**NOTE:** If you have aligner (like BWA-mem) that produses both soft- and hard-clipped reads you should use option to use only soft-clipped reads or multiply frequency assessments by 2.

### Run iJump

