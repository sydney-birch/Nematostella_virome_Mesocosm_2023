#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=bowtie
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00

##Load modules
module load bowtie2/2.5.1

## Running slurm in 3rd ali dir

## First build index of genome - built index previously
bowtie2-build -f GCF_932526225.1_jaNemVect1.1_genomic.fna Nematostella_genome

## Next run script that runs bowtie in each dir
./3.D_bowtie_ali_initial.py -b 3 -c 4th_ali_BT
