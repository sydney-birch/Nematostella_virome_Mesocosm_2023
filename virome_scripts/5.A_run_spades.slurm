#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=spades
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=16
#SBATCH --time=120:00:00

##Load modules and activate mercat2
module load spades/3.15.5

## Run spades on the filtered unmpped reads 
./5.A_run_spades.py -c ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/
