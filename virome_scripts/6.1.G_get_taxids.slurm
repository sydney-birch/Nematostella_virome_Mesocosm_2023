#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=T96_i_K_taxids
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=336:00:00

##Load modules 


## Run T0
##./4_get_taxids.sh RefSeq-release230.catalog	

## Run T96 (adjusted bash script)
##./4_get_taxids.sh RefSeq-release230.catalog

## Run T0 Clonal Individual
##./4_get_taxids.sh RefSeq-release230.catalog

## Run T96 Clonal Individual
##./4_get_taxids.sh RefSeq-release230.catalog

## Run T0 KOs Individual
##./4_get_taxids.sh RefSeq-release230.catalog

## Run T96 KOs Individual
./4_get_taxids.sh RefSeq-release230.catalog
