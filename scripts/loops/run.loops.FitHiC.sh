# #!/bin/bash
set -euo pipefail

###################################################
# Loop Calling params
###################################################
FITHIC_UPPERBOUNDS=(1000000 5000000) # default reccomendations forhuman
FITHIC_COTACTTYPES=(['cis']='intraOnly' ['trans']='interOnly' ['both']='All')
FITHIC_COTACTTYPES=(['cis']='intraOnly')
declare -rA COOLTOOLS_WEIGHTS=(["ICE"]="weight" ["Raw"]="")

###################################################
# Functions
###################################################
help() {
    echo "USAGE: $(basename ${0}) [OPTIONS] {METHOD} sample1.mcool sample{2..N}.mcool
        -a | --anaconda-dir        # where conda is installed
        -r | --resolution          # resolutions at which to annotate TADs
        -o | --output-dir          # root results dir
        -t | --ntasks-per-node
        -h | --help                # print this message
"
    exit 0
}

activate_conda() {
    # activate conda env with specific tools for each task
    case "$1" in
        fithic)    env_name="cooltools" ;;
        cooltools) env_name="cooltools" ;;
        *)         echo "Invalid conda env: $1" && exit 1 ;;
    esac
    echo "source ${CONDA_DIR}/etc/profile.d/conda.sh; conda activate ${env_name}"
}

get_sample_ID() {
    sample_ID="$(basename "${sample_file}")"
    sample_ID="${sample_ID%%.mcool}"
    sample_ID="${sample_ID%%.hic}"
    echo "${sample_ID}"
}

run_fithic() {
    sample_file="${1}"
    sample_ID="$(get_sample_ID "${sample_file}")"
    fragments_file="${2}"
    bias_file="${3}"
    resolution="${4}"
    contact_type="${5}"
    upperbound="${6}"
    lowerbound="$(echo "2 * ${resolution}" | bc)" # reccomended default
    fithic \
        --passes 1                      \
        --nOfBins 100                   \
        --biasLowerBound 0.5            \
        --biasUpperBound 2              \
        --mappabilityThres 1            \
        --outdir "${OUTPUT_DIR}"        \
        --lib "${sample_ID}"            \
        --contactType "${contact_type}" \
        --lowerbound "${lowerbound}"    \
        --upperbound "${upperbound}"    \
        --resolution "${resolution}"    \
        --biases "${bias_file}"         \
        --fragments "${fragments_file}" \
        --interactions "${sample_file}"
}

run_cooltools() {
    sample_file="${1}"
    sample_ID="$(get_sample_ID "${sample_file}")"
    resolution="${2}"
    cooltools dots "${sample_file}"
}

main() {
    # Activate conda env now if running locally (not submitting jobs)
    [[ ${EVAL_METHOD} == 'inplace' ]] && eval "${CONDA_ENV_CMD}"
    # Loop over all samples at all resolutions and call tads
    for resolution in ${RESOLUTIONS[@]}; do
        for sample_file in ${HIC_SAMPLES[@]}; do
            sample_ID="$(get_sample_ID "${sample_file}")"
            echo "${sample_ID}"
            case ${METHOD} in
                fithic)    run_fithic "${sample_file}" "${resolution}" ;;
                cooltools) run_cooltools "${sample_file}" "${resolution}" ;;
                *)         echo "Invalid method: ${METHOD}" && exit 1 ;;
            esac
        done
    done
}

###################################################
# Handle CLI args
###################################################
# Default script arguments
# BASE_DIR="/data/talkowski/Samples/16p_HiC"
BASE_DIR="./"
OUTPUT_DIR="${BASE_DIR}/results/Loops"
# EVAL_METHOD='txt'
RESOLUTIONS=(100000 50000 25000 10000)
THREADS=2
CONDA_DIR="${HOME}/miniforge3"
[[ $# -eq 0 ]] && echo "No Args" && exit 1
while getopts "ho:l:r:n:c:q:m:a:t:e:" flag; do
    case ${flag} in 
        a) CONDA_DIR="${OPTARG}" ;;
        r) RESOLUTIONS=($(echo "${OPTARG}" | cut -d',' --output-delimiter=' ' -f1-)) ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        e) EVAL_METHOD="${OPTARG}" ;;
        h) help && exit 0 ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))

###################################################
# Main 
###################################################
METHOD="${1}"
HIC_SAMPLES=${@:2}
echo "
Using Loop caller:      ${METHOD}
Using resolution(s):    ${RESOLUTIONS[*]}
Using output directory: ${OUTPUT_DIR}
Threads:                ${THREADS}"
CONDA_ENV_CMD="$(activate_conda "${METHOD}")"
OUTPUT_DIR="$(readlink -e "${OUTPUT_DIR}")"
mkdir -p "${OUTPUT_DIR}"
main 

