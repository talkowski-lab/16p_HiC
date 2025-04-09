#!/bin/bash
set -euo pipefail

BASE_DIR="/data/talkowski/Samples/16p_HiC/qc.test"
OUTPUT_DIR="${BASE_DIR}/fanc_pairs"
GENOME="${BASE_DIR}/Hi-c_ref/hg38.fasta"
CHROM_SIZES="${BASE_DIR}/Hi-c_ref/hg38.reduced.chrom.sizes"
OUTPUT_DIR="${1}"
BAM_FILES=${@:2}
for sample_file in ${BAM_FILES[@]}; do
    sample_file="$(readlink -e ${sample_file})"
    sample_ID="$(basename "$sample_file")"
    sample_ID="${sample_ID%%.lane1.hg38.0.bam}"
    pairs_file="${OUTPUT_DIR}/${sample_ID}.pairs.gz"
    [[ -e ${pairs_file} ]] && gzip -d "${pairs_file}"
    # Produce unfiltered pairs file 
    echo $sample_ID
    # pairtools parse -c ${CHROM_SIZES} ${sample_file} |
    #     pairtools sort --nproc 9 -o "${pairs_file}" 
    # fanc pairs \
    #     -f \
    #     -g "${GENOME}" \
    #     --batch-size 10000 \
    #     -r DpnII,HinfI \
    #     -t 2 \
    #     -s "${OUTPUT_DIR}/${sample_ID}.stats" \
    #     --statistics-plot "${OUTPUT_DIR}/${sample_ID}.stats.pdf" \
    #     --re-dist-plot "${OUTPUT_DIR}/${sample_ID}.re-dist.pdf" \
    #     --ligation-error-plot "${OUTPUT_DIR}/${sample_ID}.ligation-error.pdf" \
    #     ${pairs_file%%.gz} \
    #     "${pairs_file%%.gz}.fanc"
    fanc pairs \
        -f \
        --batch-size 100000 \
        -g "${GENOME}" \
        -r DpnII,HinfI \
        -t 2 \
        -s "${OUTPUT_DIR}/${sample_ID}.stats" \
        "${pairs_file%%.gz}.fanc" \
        "${pairs_file%%.gz}.filtered.fanc"
done
