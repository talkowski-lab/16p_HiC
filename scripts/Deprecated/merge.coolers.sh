#!/bin/bash

# Usage
# merge.coolers.sh 16pDELNSCHIC_S1S2.hg38.mapq_30.1000 16pDEL*.mapq_30.1000.cool
# merge.coolers.sh 16pDUPNSCHIC_S3S4.hg38.mapq_30.1000 16pDUP*.mapq_30.1000.cool
# merge.coolers.sh 16pWTNSCHIC_S5S6.hg38.mapq_30.1000 16pWT*.mapq_30.1000.cool


DISTILLER_DIR="/data/talkowski/broadIncoming/22LCC2LT4/fastq"
OUTPUT_DIR="${DISTILLER_DIR}/results/coolers_library"
RESOLUTIONS='10000000,5000000,2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000'
# Args are cool files to dump
MERGED_FILE="${1}"
HIC_SAMPLES=${@:2}
[[ -f "$MERGED_FILE" ]] && continue
cooler merge "${MERGED_FILE}.cool" ${HIC_SAMPLES[@]}
cooler zoomify \
    --out "${MERGED_FILE}.mcool" \
    --resolutions "${RESOLUTIONS}" \
    ${balance_flag} \
    "${MERGED_FILE}.cool"

