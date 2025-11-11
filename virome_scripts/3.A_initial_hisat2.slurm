#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=hisat2_intial
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00

##Load modules
module load hisat2
module load hdf5
module load sra-tools/2.10.5

## Running slurm in 1st_ali dir

## First build index of genome - built index previously - using same index as first run in: Nematostella_transcriptomics/3_HISAT2/1st_ali_11-9-23
## Run this slurm and script in Nematostella_transcriptomics/3_HISAT2 - Make a new 1st ali dir for next round of sequencing - run script there 
hisat2-build GCF_932526225.1_jaNemVect1.1_genomic.fna Nematostella_genome

## Align reads to genome - Paired reads - unzipped files
./3.A_initial_hisat2_v2.py -a ../../2_trimmomatic -b 1 -c 2nd_ali -d ../3_HISAT2/1st_ali_1-23-25/
