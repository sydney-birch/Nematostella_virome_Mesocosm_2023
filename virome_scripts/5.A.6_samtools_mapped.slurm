#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=samtools_mapped
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00

##Load modules
module load samtools

## Run sam tools
./5.A.6_samtools_mapped.py -a T0_total_mapping_fastqs 
