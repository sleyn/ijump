# iJump
Software for search of rearrangements in population sequencing data

**iJump** searches for IS elements rearrangements in evolved populations of single organism and estimates what fraction of a population is affected by the rearrangement.

## Motivation

Genome rearrangements are powerful powerfull tools for evolution in all domains of life. With current rapid decrease of sequencing costs and introduction of experimental evolution settings into laboratories more population sequencing data will be generated. It is extremely helpful for the research to see not only SNPs and short InDels but also other large rearrangements. Currently very few software tools exists that can work with data of mixed populations, not the clonal ones. The only tool for bacterial population data that was found is a [*breseq*](http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing). However although it is an excellent and comprehancise tool we've found it very slow on high coverage data.

In our experiments we found that often rearrangements are happend by mobile elements. To speed up the program we focused specifically on IS mobile elements.

## Installation

**iJump** do not need special installation but it is dependent on several Python libraries:
* **biopython**

    [Installation guidelines](https://biopython.org/wiki/Download)

* **pandas**

    [Installation guidelines](https://pandas.pydata.org/pandas-docs/stable/install.html)
    
* **pysam**

    [Installation guidelines](https://pysam.readthedocs.io/en/latest/installation.html)
    
 ## Usage
 
 ### Input
 
 iJump requires four files for input:
 1. File with mobile elements coordinates
 2. Reference DNA contigs fasta file.
 3. GFF file with reference genome annotations.
 4. BAM file of aligned Illumina reads.
 
 #### Mobile elements coordinates file
 
 File with mobile elements coordinates shoud be tab-separated tables of the following structure:
 ```
 ME_Name    Contig_Name Start_position  End_position
 ```
 
 For example:
 ```
 ISAcsp3	NODE_1	2980551	2981283
 ```

If you don't have file with coordinates of mobile elements you can find them from ISFinder website using their [BLAST](https://isfinder.biotoul.fr/blast.php) against your reference contigs.

It will return you  html page of hits that you can download and parse with **isfinder_parser.py**:
```
isfinder_parse.py -i <ISfinder BLAST HTML page>
```

### Run iJump

iJump run with the following command:
 ```
 python3 isjump.py -f <BAM file> -r <Reference FNA file> -g <Reference GFF file> -i <IS coordinates>
 ```
