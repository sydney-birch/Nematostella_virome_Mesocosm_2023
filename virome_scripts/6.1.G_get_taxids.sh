#! /bin/bash

# This bash script will run a python script fo reach file in the given dir

refseq_db=$1

for i in ./T96_KOs_50_pi_hit1_ALL_accessions_taxids/*hits
do
  echo "Getting TaxIDs for:  $i ..."
  ./4_get_taxids.py -i $i -d $1 -b $i

done

#only run mkdir once
#mkdir accid_taxid_files
#mkdir taxid_files

#mv ./T96_50_pi_hit1_SORTED_accessions/*-taxids_accids accid_taxid_files
#mv ./T96_50_pi_hit1_SORTED_accessions/*-taxids taxid_files
