#! /bin/bash

for i in ./blastout/*blastout
do
  echo "trimming blastout to 50%  $i ..."
  ./1.B_trim_blastout_50_pident.py -i $i

done

mkdir 50_pi_blastout

mv ./blastout/*50_pi ./50_pi_blastout
