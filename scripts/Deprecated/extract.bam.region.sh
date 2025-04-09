#!/bin/bash
set -euo pipefail

BASE_DIR="/data/talkowski/Samples/16p_HiC"
DELETIONS_BED="${BASE_DIR}/references.files/HiC_Deletion_Regions+1MB.bed"
# OUTPUT_DIR="${BASE_DIR}/results/Deletion.Region.BAMs"
OUTPUT_DIR="${1}"
BAM_FILES=${@:2} 

module load samtools
module load bedtools

for bam_file in ${BAM_FILES[@]}; do
    bam_file="$(readlink -e "${bam_file}")"
    sample_ID="$(basename "$bam_file")"
    sample_ID="${sample_ID%%.bam}"
    output_bam="${OUTPUT_DIR}/${sample_ID}.HiC_Deletion_Regions+1MB.bam"
    echo $sample_ID 
    time intersectBed -abam "${bam_file}" -b "${DELETIONS_BED}" >| "${output_bam}"
    time samtools sort "${output_bam}" -o "${output_bam}"
    time samtools index "${output_bam}"
done
