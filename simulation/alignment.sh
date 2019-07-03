#!/usr/bin/bash 

for i in Fastq/EC*.fq1.fq
do
    j=${i/.fq1.fq/}
    bwa_sam="$j.sam"
    bwa_bam="$j.bam"
    r2=${i/fq1/fq2}
    rgid=$(echo $j)
    rgid=${rgid/_1_trimmed.fq.gz/}
    rg="@RG\\tID:$rgid\\tSM:S1\\tPL:illumina\\tLB:L001\\tPU:R1"

    echo $rg
    ~/ngsbin/bwa-0.7.13/bwa mem -Y -t 14 -M -R "$rg" ~/Tasks/01_GEAR/Genomes/Escherichia_coli_BW25113.fna $i $r2 > $bwa_sam
    samtools sort -o $bwa_bam $bwa_sam
    samtools index -b $bwa_bam
    rm $bwa_sam
done
                                                                                                                                                                                                                                            
samtools merge Sample.bam ./Fastq/EC*.bam
samtools sort -o Sample_i.bam Sample.bam
rm Sample.bam
samtools index -b Sample_i.bam

