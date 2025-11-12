#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=samtools
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules
module load samtools

## Run sam tools - number in sam alignement file to use in samtools(if using second alignment sam file - use 2)
./4.D_samtools_processing.py -a rRNA_filtered_unmapped_fastq_files
