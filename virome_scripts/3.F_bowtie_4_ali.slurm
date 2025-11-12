#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=bowtie
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules
module load bowtie2/2.5.1

## Running slurm in 5th ali dir

## First build index of genome - built index previously
##bowtie2-build -f GCF_932526225.1_jaNemVect1.1_genomic.fna Nematostella_genome

## Next run script that runs bowtie in each dir
./3.F_bowtie_next_ali.py -b 4 -c 5th_ali_BT
