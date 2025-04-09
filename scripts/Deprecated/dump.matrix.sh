#!/bin/bash
#
DISTILLER_DIR="/data/talkowski/broadIncoming/22LCC2LT4/fastq"
OUTPUT_DIR="${DISTILLER_DIR}/results/sparse.matrices"
RESOLUTIONS=(100000 50000 500000 10000 5000)
REGIONS=(
    "chr1"
    "chr2"
    "chr3"
    "chr4"
    "chr5"
    "chr6"
    "chr7"
    "chr8"
    "chr9"
    "chr10"
    "chr11"
    "chr12"
    "chr13"
    "chr14"
    "chr15"
    "chr16"
    "chr17"
    "chr18"
    "chr19"
    "chr20"
    "chr21"
    "chr22"
    "chrX"
)
# Args are mcool files to dump
HIC_SAMPLES=${@}
for sample_file in ${HIC_SAMPLES[@]}; do
    [[ -a $sample_file ]] || continue
    sample_file="$(readlink -e "$sample_file")"
    sample_ID="$(basename "$sample_file")"
    sample_ID="${sample_ID%%.mcool}"
    for resolution in ${RESOLUTIONS[@]}; do
    for region in ${REGIONS[@]}; do
        output_file="${OUTPUT_DIR}/${sample_ID}.${resolution}.${region}.txt"
        [[ -f "$output_file" ]] && continue
        basename "${output_file}"
        cooler dump \
            --range "${region}" \
            -t pixels \
            --join \
            --no-balance \
            --out "${output_file}" \
            "${sample_file}::resolutions/${resolution}"
    done
    done
done
