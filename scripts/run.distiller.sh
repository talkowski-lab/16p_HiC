#!/bin/bash
set -euo pipefail
# Locations
DISTILLER_DIR="/data/talkowski/Samples/16p_HiC/distiller-nf"
DISTILLER_FILE="${DISTILLER_DIR}/distiller.nf"
BASE_DIR="$(pwd)"
CONDA_ENV="dist2"

# Functions
help() {
    echo TODO
}
main() {
    mkdir -p "${LOG_DIR}"
    for yml_file in ${@}; do
        sample_ID="$(basename "${yml_file}")"
        sample_ID="${sample_ID%%.distiller.yml}"
        echo "${sample_ID}"
        cmd="bash -l -c module load wget; module unload java; module load singularity/3.7.0; source "${CONDA_DIR}/etc/profile.d/conda.sh"; ${CONDA_DIR}/envs/${CONDA_ENV}/bin/nextflow run ${DISTILLER_FILE} -params-file ${yml_file} -w ${WORK_DIR} -c ${CONFIG_DIR}"
        echo "${cmd}"
        # continue 
        sbatch                                                 \
            --partition "${PARTITION}"                         \
            --time "${TIME}"                                   \
            --mem "${MEM}"                                     \
            --cpus-per-task "${CPUS}"                          \
            --ntasks-per-node "${NTASKS_PER_NODE}"             \
            --job-name "${sample_ID}-distiller"                \
            --output   "${LOG_DIR}/${sample_ID}-distiller.out" \
            --error    "${LOG_DIR}/${sample_ID}-distiller.err" \
            --wrap="${cmd}"
    done
}
# Default args
CONFIG_DIR="${DISTILLER_DIR}/nextflow.config"
WORK_DIR="${BASE_DIR}/work"
CONDA_DIR="$(conda info --base)"
LOG_DIR="${BASE_DIR}/slurm.logs"
MEM="40G"
NTASKS_PER_NODE=8
CPUS=8
PARTITION="bigmem"
TIME="4-23:59:59"
# Handle CLI args
[[ $# -eq 0 ]] && echo "No Args" && exit 1
while getopts "d:w:a:l:m:n:c:p:h" flag; do
    case ${flag} in 
        d) CONFIG_DIR="${OPTARG}" ;;
        w) WORK_DIR="${OPTARG}" ;;
        a) CONDA_DIR="${OPTARG}" ;;
        l) LOG_DIR="${OPTARG}" ;;
        m) MEM="${OPTARG}G" ;;
        n) NTASKS_PER_NODE="${OPTARG}" ;;
        c) CPUS="${OPTARG}" ;;
        p) PARTITION="${OPTARG}" ;;
        h) help && exit 0 ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))
case $PARTITION in
    short)  TIME="2:59:59"     ;;
    normal) TIME="4-23:59:59"  ;;
    long)   TIME="15-23:59:00" ;;
    # bigmem) TIME="7-23:59:59"  ;;
    bigmem) TIME="4-1:00:00"  ;;
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
# module unload java; module load wget singularity/3.7.0 source "${CONDA_DIR}/etc/profile.d/conda.sh"; conda activate dist2
main ${@}
# $(conda info --base)/envs/dist2/bin/nextflow run /data/talkowski/Samples/16p_HiC/distiller-nf -w work -params-file
