#!/bin/bash

j=$1
vars=1
for i in `seq 1 1 5`
do
    echo "bowtie ${j}" $(bowtie/bowtie-align-l --norc -k 10 -v 0 -t  indexes/chr19.${vars}.fa reads/chr19.${vars}.fa.seeds_${j}.fastq 2>&1 > /dev/null | grep time | cut -d':' -f 2,3,4 | awk 'BEGIN {FS=":"} {print $1*60*60 + $2*60 + $3}')
done
