#!/bin/bash
REF_FASTA=$1
VCF_FILE=$2
SAMPLES=$3
REGION=$4

cut -f 1 $SAMPLES | while read sample;
do 
    for i in `seq 1 1 2`
    do
        if [ -z $REGION ] 
        then
            echo ">$sample.$i"
            bcftools consensus -f $REF_FASTA -s $sample -H $i $VCF_FILE | tail -n +2
        else
            echo ">$sample.$i"
            samtools faidx $REF_FASTA $REGION | bcftools consensus -s $sample -H $i $VCF_FILE | tail -n +2
        fi
    done
done
