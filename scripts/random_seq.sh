#!/bin/bash

# usage: random_seq.sh NUM_CHARS SEQ_NAME
echo ">${2:-seq}"
echo $(cat /dev/urandom | tr -dc "ACGT" | head -c ${1:-100})
