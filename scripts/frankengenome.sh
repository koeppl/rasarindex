#!/bin/bash
INDEX=${1} # prefix of index file
READS=${2} # fq file
OUTPRE="${3:-.}"
NAME="${4:-$(basename $OUTPRE)}"

# PARAMETERS
MAX_RANGE=200
HAP=".1"
DEF_REF="hg19"
READ_LEN=100
SEQ_LEN=6000001 
FRAC=1.0

# SCRIPT LOCATIONS
RI_ALIGN=ri-align
CONVERT_OFFSETS=~/spring2018/rlfm/r-index/scripts/convert_offsets.py
FILTER_SAM=~/spring2018/rlfm/r-index/scripts/filter_sam.py
FRANKENGENOME=~/spring2018/rlfm/r-index/scripts/frankengenome.py
FKGM_TO_FA=~/spring2018/rlfm/r-index/scripts/frankengenome_to_fasta.sh

## align the reads and convert offsets

# TODO: take a subset of $READS instead
ALNS=${OUTPRE}.rlfm.sam
$RI_ALIGN --max-range $MAX_RANGE pw_locate $INDEX $READS \
    | python $CONVERT_OFFSETS --sam - ${INDEX}.fai \
    | python $FILTER_SAM --min_aln_length 21 - > $ALNS 



### create frankengenome

FKGM=${OUTPRE}.fkgm
python $FRANKENGENOME -t $HAP -d $DEF_REF -l $READ_LEN -N 6000001 $ALNS > $FKGM
$FKGM_TO_FA $INDEX $FKGM $NAME > $FKGM.fa

### build bowtie index

bowtie2-build $FKGM.fa $FKGM.fa
