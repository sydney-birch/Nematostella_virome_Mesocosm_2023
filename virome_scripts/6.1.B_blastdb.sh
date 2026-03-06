#!/bin/bash
module load blast

for i in ./fastas_prot/*.fa
do
    echo "formating for BLAST  $i ..."
    makeblastdb -in $i -parse_seqids -dbtype prot
done

# move blastdbs to new dir
mkdir blastdb
cp ./fastas_prot/* blastdb

# clean up fastas dir
rm ./fastas_prot/*.fa.*

#db types: nucl or prot
