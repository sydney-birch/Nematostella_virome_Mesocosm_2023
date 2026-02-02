#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=salmon
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=126:00:00

##Load modules 
module load salmon

## Make index of uk transcriptome
salmon index -t uk_transcriptome.fa -i mapped_index -k 31 

## Run salmon on the mapped Nematostella reads
./8_run_salmon.py -a ../4_RNA_filt/mapped_fastq_files -b ../../8_salmon -c aligned
