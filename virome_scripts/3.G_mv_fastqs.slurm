#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=fastqs
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules


## Run sam tools - number in sam alignement file to use in samtools(if using second alignment sam file - use 2)
./3.G_mv_fastqs_v2.py -b 4 -c 4_RNA_filt
