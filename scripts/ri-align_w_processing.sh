#!/bin/bash

module load python/3.6.0

PREFIX=$1 # make sure that a corresponding faidx file exists
READS=$2
MAX_HITS=$(( $3 + 1 ))
RANGE_THRES=$4

MY_DIR=$(readlink -f "$(dirname "$0")")
${MY_DIR}/../build/ri-align --max-hits $MAX_HITS --range-thres $RANGE_THRES $PREFIX.ri $READS | \
    python3 ${MY_DIR}/convert_offsets.py --sam - ${PREFIX}.fai
    # | python3 ${MY_DIR}/sam_remove_repeats.py -
