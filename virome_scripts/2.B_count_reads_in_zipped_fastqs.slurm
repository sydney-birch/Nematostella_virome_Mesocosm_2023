#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=count_reads
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules

## Run count reads script 
## trimmed reads after trimmomatic
./2.B_count_reads_in_zipped_fastqs.py -a 2_trimmomatic

## Raw reads 5-17-24 
##./count_reads_in_zipped_fastqs.py -a raw_reads_5-17-24

## Raw reads 10-7-24
##./count_reads_in_zipped_fastqs.py -a raw_reads_10-7-24
