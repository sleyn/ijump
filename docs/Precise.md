# "Precise" workflow

With the "Precise" workflow iJump tries to find the precise location of IS element insertion.

## Algorithm

The algorithm is following:

1. Make a forward search of insertions:
    a. Collect clipped reads near the boundaries of IS elements using provided IS elements coordinates.
    b. BLAST clipped parts of the reads against reference genome (parts of the reads that are not aligned to the reference).
    c. Collect positions of junctions from hits if:
        - Hit is a single best (filter out repeat regions).
        - Hit is NOT NEAR any IS elements boundaries.
2. Make a reverse search:
    a. Collect clipped reads near found positions of junctions.
    b. BLAST clipped parts of the reads against reference genome.
    c. Collect positions of junctions from hits if:
        - Hit is a single best (filter out repeat regions).
        - Hit is NEAR any IS elements boundaries.
    d. Using this information assign IS elements back to the reads. 
3. For each IS element if it is possible find a pair of junctions positions to the left and right part of the IS element that have:
    a. Minimal difference in read counts.
    b. Have distance less than 30bp.
    This positions will go to the output.
4. For junction position take back number of reads supporting it (from the step 2a) and recalculate number of reads that are assigned to IS element using proportional counts.
    For example:
    - 500 reads were found for position 1234.
    - By reverse search 300 reads were assigned to *IS1* and 100 reads to *IS2*.
    - Using proportions: $500 * 300 / (300 + 100) = 375$ reads will be assigned to *IS1* and $500 * 100 / (300 + 100) = 125$ reads will be assigned to *IS2*.
5. iJump calculates frequency of insertion using average of frequencies of two junctions. For each junction frequency estimation we use several numbers:
    - $N_{cl}$ - Number of clipped reads that support junction position for a given IS element.
    - $N_{ov}$ - Number of clipped reads that overlap with junction but are clipped at different coordinate.
    - $N_{ncl}$ - Number of non-clipped reads that overlap with the junction position.
    - $C_{match}$ - Correction coefficient that accounts to inability of aligner to map short clipped portions of reads. This causes loss of observed coverage.
        $C_{match} = \frac{1}{1 - m_{match} / a_{rlen}}$
        where:
            - $a_{rlen}$ - average read length
            - $m_{match}$ - minimum length of the read part that could be aligned to reference. Accessed as the minimum of longest clipped part of the read (*e.g.* for read with CIGAR string 10S120M30S *$m_{match}$* is 30).
      
    The frequency will be calculated by formula:
    $ Freq = \frac{N_{cl} * C_{match}}{N_{cl} * C_{match} + N_{ov} * C_{match} + N_{ncl} + 0.1}

### Output

#### ijump_junction_pairs.txt

The main results output. File contains information about insertions. File contains following columns:

* *Position_[l|r]*<br>
	junction coordinate

* *Chrom*<br>
    	contig name where junction was found

* *Count_mapped_to_IS_[l|r]*<br>
    	number of reads that support junction to the IS elements

* *IS_name*<br>
    name of the IS element
  
* *N_clipped_[l|r]*<br>
    number of clipped reads supporting junction
  
* *N_clipped_[l|r]_corrected*<br>
    number of clipped reads corrected by loss of reads that were not mapped due to the small clipped part
  
* *N_clipped_[l|r]_correction*<br>
    estimated number of clipped reads not mapped to the junction coordinate
  
* *Depth*<br>
    average depth at junction coordinates<br>
    If only one coordinate is present - the value is the depth at this coordinate.

* *Dist*<br>
    distance between junctions

* *N_overlap_formula_[l|r]*<br>
    a number of clipped reads that are ovelapped with junction position but have another junction coordinate used for frequency estimation.<br>
    Different from *N_overlap_[l|r]* as it is equal *N_overlap_[l|r]_corrected - N_clipped_[l|r]_corrected* *Position_r > x.Position_l* and *N_overlap_[l|r]_corrected* otherwise.

* *N_overlap_[l|r]*<br>
    a number of clipped reads that are ovelapped with junction position but have another junction coordinate

* *N_overlap_[l|r]_corrected*<br>
    *N_overlap_[l|r]* corrected by loss of reads that were not mapped due to the small clipped part

* *N_overlap_[l|r]_correction*<br>
    estimated number of clipped reads not mapped to the region that would overlap with junction

* *N_unclipped_[l|r]*<br>
    number of not clipped reads that overlap with junction 
  
* *Frequency*<br>
    estimated frequency of insertion
  
* *Frequency_[l|r]*<br>
    estimated frequency of insertion for each junction


#### ijump_junctions.txt

Contains information about junctions for each read. File contains following columns:

* *index*  
	 order number

* *ID*  
	 unique identifier

* *IS name*  
	 mobile element name

* *IS pos*  
	 what part of the read matches mobile element

* *IS chrom*
	 name of contig where mobile element is located in the reference

* *Read name*
	 read name where junction was observed

* *Chrom*
	 name of contig where mobile element jumped

* *Position*
	 posistion of the junction

* *Orientation*
	 orientation of mobile element relative to junction

* *Note* 
	 mark if junction is in other mobile elements - usually indicates false positive hits

* *Locus tag*
	 locus tag of the affected gene; in the case of intergenic region two locus tags will be shown with us_ or ds_ prefixes that indicate upstream or downstream position of the region relative to the genes.

* *Gene*  
	 trivial name of the affected gene
  
