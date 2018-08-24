#!/bin/bash

REF=$1
FKGM=$2
NAME=$3
module load samtools

echo ">$NAME"
cat $FKGM | xargs samtools faidx $1 | grep -v ">" | tr -d '\n' | fold -w 60 
