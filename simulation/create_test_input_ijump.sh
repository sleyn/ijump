#!/bin/bash

#set up enviroment for ijump simulation
conda create -y -n ijump_simulation
conda activate ijump_simulation
conda config --add channels bioconda
conda config --add channels conda-forge

#Install packages
conda install -y -c bioconda samtools
conda install -y -c bioconda bwa
apt-get install -y art-nextgen-simulation-tools

#Prepare output dir and give samtools the right to edit files
mkdir Fastq
chmod -R 777 *

#Create test input for ijump:
./makefq.sh
./alignment.sh