#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=counts
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00

##Load modules 

##Get blastout hit counts (wc -l)
##./1.A_get_blast_hit_counts.py -a blastout

## Get 70 pident counts
##./1.A_get_blast_hit_counts.py -a 70_pi_blastout

## Get counts on unique 50 pid accids
./1.A_get_blast_hit_counts.py -a 50_pi_hit1_SORTED_accessions
