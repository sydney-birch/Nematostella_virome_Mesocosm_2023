#! /bin/bash

# This bash script will run a python script fo reach file in the given dir
#$1 = input dir to use

for i in $1/*-taxids
do
  echo "Running Taxon Kit for:  $i ..."
  ./5_run_taxonkit.py -i $i 

done


mv *linage.txt ./lineage_files
