NUM_SAMPLES=${1}
ALL_VAR_FASTA=chr19_vars.fa
PREFIX=chr19.${NUM_SAMPLES}
RLFM_DIR=~/scratch/rlfm/r-index/build
BOWTIE_DIR=~/scratch/rlfm/bowtie/
INDEX_DIR=indexes
OUT_FASTA=${INDEX_DIR}/${PREFIX}.fa

# copy over the first $NUM_SAMPLES records over to a new file
awk -v NSEQS="${NUM_SAMPLES}" '/^>/ {n++} n>NSEQS {exit} {print}' ${ALL_VAR_FASTA} > ${OUT_FASTA}

samtools faidx ${OUT_FASTA}

${RLFM_DIR}/ri-buildfasta -o ${OUT_FASTA} ${OUT_FASTA}
${BOWTIE_DIR}/bowtie-build-l --threads 16 $OUT_FASTA ${OUT_FASTA}
