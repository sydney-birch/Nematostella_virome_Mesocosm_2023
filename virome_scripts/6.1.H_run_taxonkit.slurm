#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=taxonkit
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00

##Load modules 


## Run taxonkit bash script (calls a python script)
./5_run_taxonkit.sh final_taxids
