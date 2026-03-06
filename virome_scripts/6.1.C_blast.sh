#! /bin/bash
module load blast

fasta_query=$1

for fasta in ./blastdb/*.fa
do
    echo "BLASTING  $fasta ..."
    blastp -query $1 -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1
done

mkdir blastout

mv ./blastdb/*blastout ./blastout/

