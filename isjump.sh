#!/bin/bash

python3 isjump.py -f Sample_r_bqsr.bam -r A_baumannii_assembly.fna -g A_baumannii_assembly.gff
circos -config ./data/circos.conf