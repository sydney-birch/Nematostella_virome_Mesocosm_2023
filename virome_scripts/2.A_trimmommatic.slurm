#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00

##Load module
module purge
module load trimmomatic/0.39 

## Running slurm in 2_trimmomatic dir - run for each set of reads 
##./2.A_trimmomatic.py -a ../raw_reads_10-7-24 -b ../2_trimmomatic
./2.A_trimmomatic.py -a ../raw_reads_5-17-24 -b ../2_trimmomatic
	##Re-running it stopped partially through 5-17 reads - will see if picks back up - it might just restart and overright  
