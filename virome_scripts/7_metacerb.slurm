#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=metacerberus
#SBATCH --nodes=1
##SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=120:00:00

##Load modules and activate mercat2
module load anaconda3
conda activate metacerberus

## Official Run - individual nuc assemblies
metacerberus.py --prodigal ../5_assemblies/fastas_nuc  --hmm ALL --dir_out metacerb_indiv_nuc_2-21-25

rm -r tmp-hydra/

## Official Run - collapsed nuc assemblies
metacerberus.py --prodigal ../5_assemblies/fastas_nuc_collapsed  --hmm ALL --dir_out metacerb_collapsed_nuc_2-21-25




## TEST: combine reps in one run from individual nuc assembly - worked
##metacerberus.py --prodigal SC_T0-WT_test --hmm ALL --dir_out indiv_T0-WT_Test

## TEST: compare to collapsed rep -worked 
##metacerberus.py --prodigal ../5_assemblies/fastas_nuc_collapsed/SC_T0-WT-nuc.fa  --hmm ALL --dir_out coll_T0-WT_Test

## TEST: combine reps from prot individual assemblies .fa extensions - doesn't work - registers files but couldn't read seqs 
##metacerberus.py --prodigal test_SC_T0-WT_prot --hmm ALL --dir_out prot_test_1_indiv_T0-WT

## TEST: combine reps from prot individual assemblies .faa extensions - doesn't work - doesn't register files
##metacerberus.py --prodigal test_SC_T0-WT_prot --hmm ALL --dir_out prot_test_2_indiv_T0-WT
