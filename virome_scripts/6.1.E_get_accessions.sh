#! /bin/bash

for i in ./50_pi_blastout/*50_pi
do
  echo "getting hit accessions  $i ..."
  cut -f1 $i > ${i%.}_hits
done

mkdir 50_pi_hit1_accessions

mv ./50_pi_blastout/*hits ./50_pi_hit1_accessions
