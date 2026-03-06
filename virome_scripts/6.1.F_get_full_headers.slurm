#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=headers
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=48:00:00

##Load modules 
##module load blast

## Run blast
./3_get_full_headers.sh viral.1.protein.faa
