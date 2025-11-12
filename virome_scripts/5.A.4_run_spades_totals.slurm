#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=spades
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=16
#SBATCH --time=120:00:00

##Load modules and activate mercat2
module load spades/3.15.5

## Run spades on total T0  
spades.py --rnaviral -1 SC_T0-total_R1.fq.gz -2 SC_T0-total_R2.fq.gz -o SC_T0_Total_spades_output

## Run spades on total T96
spades.py --rnaviral -1 SC_T96-total_R1.fq.gz -2 SC_T96-total_R2.fq.gz -o SC_T96_Total_spades_output
