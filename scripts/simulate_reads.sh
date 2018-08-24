echo $1 $2 $3 $4 >> errs
REF=$1
PREFIX=$2
LEN=$3
NUM_READS=$4
MUT_RATE=${5:-0.0}
DIR=$(readlink -f "$(dirname "$0")")

# $DIR/wgsim -N $4 -1$3 -2$3 -h -e $MUT_RATE -r 0.0 -R 0.0 $1 $2.fq /dev/null
mason illumina -snN -nN -N ${NUM_READS} -f -sq -n ${LEN} -pi 0 -pd 0 --pmm $MUT_RATE -pmmb $MUT_RATE --pmme $MUT_RATE  -hs 0 -hi 0 -hnN -o ${PREFIX}.fq ${REF} 

