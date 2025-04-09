#!/bin/bash
# set -euo pipefail

BASE_DIR="/data/talkowski/Samples/16p_HiC"
DISTILLER_FILE="${BASE_DIR}/distiller-nf/distiller.nf"
LOG_DIR="${BASE_DIR}/slurm.logs"
YML_FILE="$1"
RUN_NAME="$(basename "$1")"
RUN_NAME="${RUN_NAME%%.distiller.yml}"

sbatch \
    --partition bigmem \
    --nodes=1 \
    --ntasks=8 \
    --cpus-per-task=8 \
    --mem=40G \
    --time 3-1:00:00 \
    --job-name "${RUN_NAME}-distiller" \
    --output   "${LOG_DIR}/${RUN_NAME}-distiller.out" \
    --error    "${LOG_DIR}/${RUN_NAME}-distiller.err" \
    --wrap="bash -l -c module load wget
module unload java
module load FastQC/0.11.8-Java-1.8
module load samtools/1.11
module singularity/3.7.0
module load bwa/0.7.17
source ${HOME}/miniforge3/etc/profile.d/conda.sh
conda activate dist2
${HOME}/miniforge3/envs/distiller/bin/nextflow run ${DISTILLER_FILE} -params-file ${YML_FILE} -resume"

