#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=spades
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=16
#SBATCH --time=120:00:00

##Load modules and activate 
module load transdecoder/5.7.1

## Run transdecoder on each assembly (in their own dirs) 
./5.B.2_run_transdecoder_collapse_reps.py -c fastas_nuc_collapsed -b fastas_prot_collapsed
