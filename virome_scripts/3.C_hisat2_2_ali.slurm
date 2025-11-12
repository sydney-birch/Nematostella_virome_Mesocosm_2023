#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=hisat2_2
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00

##Load modules
module load hisat2
module load hdf5
module load sra-tools/2.10.5

##Run script: 
./3.C_hisat2_next_ali_v2.py -b 2 -c 3rd_ali_BT
