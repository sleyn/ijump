# iJump
Software for search of rearrangements in population sequencing data

**iJump** analizes BAM files of aligned Illumina short reads for population sequencing data. It searches for IS elements rearrangements and estimates what fraction of a population is affected by the rearrangement.

## Motivation

Genome rearrangements are powerful powerfull tools for evolution in all domains of life. With current rapid decrease of sequencing costs and introduction of experimental evolution settings into laboratories more population sequencing data will be generated. It is extremely helpful for the research to see not only SNPs and short InDels but also other large rearrangements. Currently very few software tools exists that can work with data of mixed populations, not the clonal ones. The only tool for bacterial population data that was found is a [*breseq*](http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing). However although it is an excellent and comprehancise tool I've found it very slow on high coverage data.

In our experiments we found that often rearrangements are happend by mobile elements. To speed up the program we focused specifically on IS mobile elements.

## Installation

**iJump** is dependent on several python libraries:
* **biopython**

    Installation guidelines could be found [https://biopython.org/wiki/Download](https://biopython.org/wiki/Download)

* **pandas**
