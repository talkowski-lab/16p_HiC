#!/bin/bash
set -euo pipefail

# /data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes
BASE_DIR="/data/talkowski/Samples/16p_HiC/qc.test"
OUTPUT_DIR="${BASE_DIR}/Hi-c_ref"
SIZES="${OUTPUT_DIR}/hg38.reduced.chrom.sizes"
REFERENCE="${OUTPUT_DIR}/hg38.fasta"
cooler digest $SIZES $REFERENCE DpnII > "${OUTPUT_DIR}/hg38.DpnII.digested.bed"
cooler digest $SIZES $REFERENCE HinfI > "${OUTPUT_DIR}/hg38.HinfI.digested.bed"
bedops --partition \
    "${OUTPUT_DIR}/hg38.DpnII.digested.bed" \
    "${OUTPUT_DIR}/hg38.HinfI.digested.bed" \
    >| "${OUTPUT_DIR}/hg38.ARIMA.digested.bed"
