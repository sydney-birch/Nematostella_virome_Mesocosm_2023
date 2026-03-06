#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00

##Load modules 
module load blast

#Run Blast
./1_blast.sh viral.1.protein.faa
