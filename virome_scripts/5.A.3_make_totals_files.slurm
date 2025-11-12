#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=total_files
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00


## Running slurm 

## First run on T0 data
./5.A.3_make_totals_files.py -a ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/SC_T0_fastqs -b filtered -c ../../../5_assemblies/ -d SC_T0

## Next run on T96 data 
./5.A.3_make_totals_files.py -a ../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/SC_T96_fastqs -b filtered -c ../../../5_assemblies/ -d SC_T96
