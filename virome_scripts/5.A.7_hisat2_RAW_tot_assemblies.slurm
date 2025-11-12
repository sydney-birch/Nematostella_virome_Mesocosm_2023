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
##./5.A.7_hisat2_RAW_tot_assemblies.py -a T0_total_RAW_mapping_fastqs -b total-SC-T0_assembly.fa_FIX.fa -c ../../raw_reads_5-17-24/ -d T0
./5.A.7_hisat2_RAW_tot_assemblies.py -a T0_total_RAW_mapping_fastqs -b total-SC-T0_assembly.fa_FIX.fa -c ../../raw_reads_10-7-24/ -d T0


# Run T96
##./5.A.7_hisat2_RAW_tot_assemblies.py -a T96_total_RAW_mapping_fastqs -b total-SC-T96_assembly.fa_FIX.fa -c ../../raw_reads_5-17-24/ -d T96
##./5.A.7_hisat2_RAW_tot_assemblies.py -a T96_total_RAW_mapping_fastqs -b total-SC-T96_assembly.fa_FIX.fa -c ../../raw_reads_10-7-24/ -d T96

