#!/bin/bash

j=$1
vars=50

# rlfm
for i in $(echo old fix_lf fix_lf2 fix_lf3 fix_lf4)
do
    for n in `seq 1 1 5`
    do
        # echo "$i $j" $(r-index/${i}/ri-time --niter 1 --max-hits 10 indexes/chr19.${vars}.fa.ri reads/chr19.${vars}.fa.seeds_${j}.pizza 2>&1 | awk '/Took/ {printf "%s ", substr($2, 1, length($2)-1)} /LF_STATS/ {printf "%s %s %s %s %s %s %s", $2, $3, $4, $5, $6, $7, $8; print"\n"}')
        echo "$i $j" $(r-index/${i}/ri-time --niter 1 --max-hits 10 indexes/chr19.${vars}.fa.ri reads/chr19.${vars}.fa.seeds_${j}.pizza 2>&1 | awk '/Took/ {printf "%s ", substr($2, 1, length($2)-1)}')
    done
done

# bowtie
for i in `seq 1 1 5`
do
    echo "bowtie ${j}" $(bowtie/bowtie-align-l --norc -k 10 -v 0 -t  indexes/chr19.${vars}.fa reads/chr19.${vars}.fa.seeds_${j}.fastq 2>&1 > /dev/null | grep time | cut -d':' -f 2,3,4 | awk 'BEGIN {FS=":"} {print $1*60*60 + $2*60 + $3}')
done
