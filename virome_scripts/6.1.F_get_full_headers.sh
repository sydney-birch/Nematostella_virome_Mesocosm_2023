#! /bin/bash

full_fasta_file=$1
 
for i in ./50_pi_hit1_SORTED_accessions/*unique_accids
do
  echo "Getting full headers for $i ...."
  ./3.B_get_full_headers.py -a $i -f $1 -o ${i%.}_total
done

mkdir 50_pi_total_header_hit1_accessions
mv ./50_pi_hit1_SORTED_accessions/*_total ./50_pi_total_header_hit1_accessions
