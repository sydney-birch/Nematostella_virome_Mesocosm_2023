#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=hisat2
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules
module load hisat2
module load hdf5
module load sra-tools/2.10.5

##Run script: 
./4.C_filtered_hisat2.py -a rRNA_filtered_unmapped_fastq_files -b rRNAdb_ed1.fa -c unmapped_fastq_files
