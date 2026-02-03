#!/bin/bash
set -uo pipefail

###################################################
# cooltools Compartment Calling params
###################################################
# declare -rA COMPARTMENT_WEIGHTS=(["ICE"]="weight" ["Raw"]="RAW")
declare -rA COMPARTMENT_WEIGHTS=(["ICE"]="weight")
REF_DIR="/data/talkowski/tools/ref/Hi_c_noalt"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
CHROMSIZES_FILE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
REFERNCE_FILES_DIR="/data/talkowski/Samples/16p_HiC/reference.files"
GENOME_BINS_DIR="${REFERNCE_FILES_DIR}/genome.bins"
TRACK_FILES_DIR="${REFERNCE_FILES_DIR}/genome.tracks"

###################################################
# Functions
###################################################
help() {
    echo "USAGE: $(basename "${0}") [OPTIONS] {MODE} sample1.mcool sample{2..N}.mcool
        -m | mode: phasing, compartment
        -a | where conda is installed
        -r | resolutions at which to annotate compartments
        -o | root results dir
        -n | number of threads
        -h | print this message"
    exit 0
}

get_sample_ID() {
    sample_ID="$(basename "${sample_file}")"
    sample_ID="${sample_ID%%.mcool}"
    sample_ID="${sample_ID%%.hic}"
    echo "${sample_ID}"
}

compute_phasing_track() {
    # Make reference files to use for calling + oritenting compartments with cooltools
    local resolution="${1}"
    # bin the genome at a specific resolution
    bins_file="${GENOME_BINS_DIR}/resolution_${resolution}/genome.bins.tsv"
    if ! [[ -e "${bins_file}" ]]; then 
        mkdir -p "$(dirname "${bins_file}")"
        echo ${bins_file}
        cooltools genome binnify "${CHROMSIZES_FILE}" "${resolution}" >| "${bins_file}"
    fi
    # calculate frac of nucleotides that are GC per bin
    track_file="${TRACK_FILES_DIR}/resolution_${resolution}/${TRACK_TYPE}-genome.track.tsv"
    if ! [[ -e "${track_file}" ]]; then 
        mkdir -p "$(dirname "${track_file}")"
        echo "${track_file}"
        cooltools genome ${TRACK_TYPE} ${bins_file} ${GENOME_REFERENCE} >| "${track_file}"
    fi
}

run_cooltools_compartments() {
    local sample_file="${1}"
    local resolution="${2}"
    sample_ID="$(get_sample_ID "${sample_file}")"
    for weight_name in "${!COMPARTMENT_WEIGHTS[@]}"; do
        for resolution in ${RESOLUTIONS[@]}; do
        # Now produce eigenvector + orient with tack file (i.e. compartment per bin)
        param_dir="${OUTPUT_DIR}/method_${METHOD}/track_${TRACK_TYPE}/weight_${weight_name}/resolution_${resolution}"
        # mkdir -p "${param_dir}"
        output_prefix="${param_dir}/${sample_ID}-"
        # List other input files i.e. list of bins + phasing track data
        bins_file="${GENOME_BINS_DIR}/resolution_${resolution}/genome.bins.tsv"
        track_file="${TRACK_FILES_DIR}/resolution_${resolution}/${TRACK_TYPE}-genome.track.tsv"
        # echo "
        # cooltools eigs-cis
        #     --phasing-trace ${track_file}
        #     --n-eigs 3
        #     --clr-weight-name ${COMPARTMENT_WEIGHTS[${weight_name}]}
        #     -o ${output_prefix}
        #     ${sample_file}
        # ====================================="
        cooltools eigs-cis --phasing-track "${track_file}" --n-eigs 3 --clr-weight-name ${COMPARTMENT_WEIGHTS[${weight_name}]} -o "${output_prefix}"                                    "${sample_file}" 
        echo "====================================="
        done
    done
}

main() {
    # Loop over all samples at all resolutions and call tads
    local hic_samples=${*}
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
    conda activate 'cooltools'
    if [[ ${MODE} == 'phasing' ]]; then
        for resolution in ${RESOLUTIONS[@]}; do
            compute_phasing_track ${resolution}
        done
    elif [[ ${MODE} == 'compartment' ]]; then
        for sample_file in ${hic_samples[@]}; do
            sample_ID="$(get_sample_ID "${sample_file}")"
            echo "${sample_ID}"
            run_cooltools_compartments "${sample_file}" ${resolution}
        done
    fi
}

###################################################
# Handle CLI args
###################################################
# Default script arguments
# BASE_DIR="/data/talkowski/Samples/16p_HiC"
BASE_DIR="./"
OUTPUT_DIR="${BASE_DIR}/results/compartments"
METHOD="cooltools"
TRACK_TYPE="genecov"
RESOLUTIONS=(100000 50000 25000 10000 5000)
CONDA_DIR="${HOME}/miniforge3"
THREADS="$(nproc)"
# Default SLURM params
# QUEUE="normal"; MEM_GB=30; NTASKS_PER_NODE=1; CPUS=2; THREADS=2; LOG_DIR="${BASE_DIR}/slurm.logs"
[[ $# -eq 0 ]] && echo "No Args" && exit 1 # while getopts "ho:l:r:n:c:q:m:a:t:e:" flag; do
while getopts "hm:a:r:o:n:" flag; do
    case ${flag} in 
        m) MODE="${OPTARG}" ;;
        t) TRACK_TYPE="${OPTARG}" ;;
        a) CONDA_DIR="${OPTARG}" ;;
        r) RESOLUTIONS=($(echo "${OPTARG}" | cut -d',' --output-delimiter=' ' -f1-)) ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        n) THREADS="${OPTARG}" ;;
        h) help && exit 0 ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))

###################################################
# Main 
###################################################
echo "
Using Phasing Track:      ${TRACK_TYPE}
Using resolution(s):      ${RESOLUTIONS[*]}
Using output directory:   ${OUTPUT_DIR}
Threads:                  ${THREADS}"
OUTPUT_DIR="$(readlink -e "${OUTPUT_DIR}")"
mkdir -p "${OUTPUT_DIR}"
main "${@}"

