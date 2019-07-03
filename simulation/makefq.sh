#!/usr/bin/bash
art_illumina \
	--seqSys HSXn \
	--len 150 \
	--mflen 550 \
	--noALN \
	--paired \
	--fcov 190 \
	--in References/EC_1.fasta \
	--sdev 10 \
	--out Fastq/EC1.fq
	
art_illumina \
	--seqSys HSXn \
	--len 150 \
	--mflen 550 \
	--noALN \
	--paired \
	--fcov 250 \
	--in References/EC_2.fasta \
	--sdev 10 \
	--out Fastq/EC2.fq

art_illumina \
	--seqSys HSXn \
	--len 150 \
	--mflen 550 \
	--noALN \
	--paired \
	--fcov 50 \
	--in References/EC_3.fasta \
	--sdev 10 \
	--out Fastq/EC3.fq
	
art_illumina \
	--seqSys HSXn \
	--len 150 \
	--mflen 550 \
	--noALN \
	--paired \
	--fcov 10 \
	--in References/EC_4.fasta \
	--sdev 10 \
	--out Fastq/EC4.fq
	
art_illumina \
	--seqSys HSXn \
	--len 150 \
	--mflen 550 \
	--noALN \
	--paired \
	--fcov 500 \
	--in References/EC_5.fasta \
	--sdev 10 \
	--out Fastq/EC5.fq
