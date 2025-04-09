#!/bin/bash
set -euo pipefail

REF_GENOME='/data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes'
THREADS=64
JUICER_JAR="/data/talkowski/Samples/WAPL_NIPBL/HiC/distiller-nf/juicer_tools_1.22.01.jar"
OUTPUT_DIR=${1}
PAIRS_FILES=${@:2}
for pairs_file in ${PAIRS_FILES[@]}; do
    sample_ID="$(basename $pairs_file)"
    sample_ID="${sample_ID%%.nodups.pairs.gz}"
    echo $pairs_file
    echo $sample_ID
    java -Xmx48000m -Djava.awt.headless=true -jar "${JUICER_JAR}" pre \
        --threads $THREADS \
        "$pairs_file" \
        "${OUTPUT_DIR}/${sample_ID}.hic" \
        "$REF_GENOME"
done
