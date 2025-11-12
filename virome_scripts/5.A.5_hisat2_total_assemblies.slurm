#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=hisat2
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules
module load hisat2
module load hdf5
module load sra-tools/2.10.5

##Run T0: 
./5.A.5_hisat2_total_assemblies.py -a T0_total_mapping_fastqs -b total-SC-T0_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T0


# Run T96
##./5.A.5_hisat2_total_assemblies.py -a T96_total_mapping_fastqs -b total-SC-T96_assembly.fa_FIX.fa -c ../../4_RNA_filt/rRNA_filtered_unmapped_fastq_files/ -d T96

