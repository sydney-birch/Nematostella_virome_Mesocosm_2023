#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=count_reads
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules

## run on nuc fastas
##./5.C_count_headers.py -a fastas_nuc
./5.C_count_headers.py -a fastas_nuc_collapsed

## Run on prot fastas
##./5.C_count_headers.py -a fastas_prot
./5.C_count_headers.py -a fastas_prot_collapsed


