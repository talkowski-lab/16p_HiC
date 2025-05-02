#!/bin/bash
set -euo pipefail
set -x
# module unload java; module load wget singularity/3.7.0 FastQC/0.11.8-Java-1.8 samtools/1.11 bwa/0.7.17; conda activate dist2
# Locations
BASE_DIR="/data/talkowski/Samples/16p_HiC"
DISTILLER_FILE="${BASE_DIR}/distiller-nf/distiller.nf"
LOG_DIR="${BASE_DIR}/slurm.logs"
mkdir -p "${LOG_DIR}"
CONDA_ENV="dist2"
#CONDA_DIR="$(conda info --base)"

CONDA_DIR="/PHShome/ry077/miniforge3"
# Functions
help() {
    echo TODO
}
main() {
    for yml_file in ${@}; do
        sample_ID="$(basename "${yml_file}")"
        sample_ID="${sample_ID%%.distiller.yml}"
        echo $sample_ID
        log_file="${LOG_DIR}/${sample_ID}"-distiller
        cmd="bash -l -c module load wget; module load samtools; module unload java; module unload perl; module load singularity/3.7.0; export PATH="/PHShome/ry077/miniforge3/:$PATH"; source ${CONDA_DIR}/etc/profile.d/conda.sh; conda activate ${CONDA_ENV}; source ${CONDA_DIR}/etc/profile.d/conda.sh; ${CONDA_DIR}/envs/${CONDA_ENV}/bin/nextflow run ${DISTILLER_FILE} -params-file ${yml_file} -w ${WORK_DIR} -c ${CONFIG_DIR}"
        sbatch \
            --partition ${PARTITION} \
            --time ${TIME} \
            --mem ${MEM} \
            --ntasks-per-node ${NTASKS_PER_NODE} \
            --cpus-per-task ${CPUS} \
            --job-name "${sample_ID}-distiller" \
            --output   "${LOG_DIR}/${sample_ID}-distiller.out" \
            --error    "${LOG_DIR}/${sample_ID}-distiller.err" \
            --mail-type "ALL" \
            --mail-user "ryadav1@mgh.harvard.edu" \
            --wrap="${cmd}"
    done
}
# Default args
THREADS=8
PARTITION="normal"
TIME="4-23:59:59"
MEM="10G"
NTASKS_PER_NODE=8
CPUS=4
WORK_DIR="${BASE_DIR}/work_RY"
CONFIG_DIR="${BASE_DIR}/distiller-nf/nextflow.config"
# Handle CLI args
[[ $? -ne 0 ]] && echo "No Args" && exit 1
VALID_ARGS=$(getopt -o ht:c:p:m:w:a: --long help,mem,partition,cpus,ntasks_per_node,work-dir,anaconda-dir -- "$@")
eval set -- "$VALID_ARGS"
while [ : ]; do
    case "$1" in
        -n|--nextflow-config)
            CONFIG_DIR="${2}"
            shift 2
            ;;
        -a|--anaconda-dir)
            CONDA_DIR="${2}"
            shift 2
            ;;
        -w|--work-dir)
            WORK_DIR="${2}"
            shift 2
            ;;
        -m|--mem)
            MEM="${2}G" 
            shift 2
            ;;
        -t|--ntasks-per-node)
            NTASKS_PER_NODE="${2}" 
            shift 2
            ;;
        -c|--cpus)
            CPUS="${2}" 
            shift 2
            ;;
        -p|--partition)
            PARTITION="${2}" 
            shift 2
            ;;
        -h|--help) 
            help 
            exit 0 
            ;;
        --)
            shift 
            break
            ;;
    esac
done
case $PARTITION in
    short)  TIME="2:59:59"     ;;
    normal) TIME="4-23:59:59"  ;;
    long)   TIME="15-23:59:00" ;;
    bigmem) TIME="7-23:59:59"  ;;
    *) echo "Invalid partiton $PARTITION" && exit 1 ;;
esac
# Print args
echo "Args to launch each job with:
PARTITION:       ${PARTITION}
TIME:            ${TIME}
MEM:             ${MEM}
NTASKS_PER_NODE: ${NTASKS_PER_NODE} 
CPUS:            ${CPUS} 
WORK_DIR:        ${WORK_DIR}"
# All args are assumed to be distiller compatible yml files
main ${@}
