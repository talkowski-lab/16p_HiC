#!/bin/bash
set -o pipefail

###################################################
# Fixed variables/lists 
###################################################
RESOLUTIONS=(100000 50000 25000 10000 5000)
declare -rA WEIGHTS=([raw]='' [balanced]='weight')

###################################################
# Functions
###################################################
help() {
    echo "Usage: ${0} [OPTIONS] ${OUTPUT_DIR} sample1.mcool sample2.mcool ...
 -o  | output_dir: path to save output files
 -t  | threads: number of threads to use
 -f  | force_redo: the script automatically detects if coverage files exists and will, by default, not recompute existing results. this flag forces re-computing even if the results exist
"
    exit 0
}

compute_matrix_coverage() {
    output_dir="$(readlink -e "${1}")"
    mkdir -p "${1}"
    echo "Saving results in ${output_dir}"
    for sample_file in "${@:2}"; do
        # Extract SampleID
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        sample_file="$(readlink -e "${sample_file}")"
        # loop over all resolutions + normalizations
        for resolution in "${RESOLUTIONS[@]}"; do
            uri="${sample_file}::resolutions/${resolution}"
        # for uri in $(cooler ls "${sample_file}"); do 
            # resolution="$(echo "${uri}" | rev | cut -d '/' -f1 | rev)"
            for weight_name in "${!WEIGHTS[@]}"; do
                weight="${WEIGHTS[${weight_name}]}"
                if [[ ${weight} != '' ]]; then
                    weight_flag="--clr_weight_name ${weight} "
                else 
                    weight_flag=""
                fi
                # name output file directory with params
                param_dir="${output_dir}/weight_${weight_name}/resolution_${resolution}"
                output_file="${param_dir}/${sample_ID}-coverage.tsv"
                # skip if output file exists and no -f flag
                if [[ -f "${output_file}" && ${FORCE_REDO} == 0 ]]; then
                    echo "Skipping, results file exists: ${output_file}"
                else
                    # Caculate total IF of each bin
                    echo "Generating results file: ${output_file}"
                    mkdir -p "${param_dir}"
                    cooltools coverage ${weight_flag}--nproc "${THREADS}" --output "${output_file}" "${uri}"
                fi
            done
        done
    done
}

###################################################
# Handle Arguments
###################################################
THREADS=$(nproc)
FORCE_REDO=0
# Handle CLI args
[[ $# -eq 0 ]] && echo "No Args" && help
while getopts "a:t:fh" flag; do
    case ${flag} in 
        o) OUTPUT_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        f) FORCE_REDO=1 ;;
        h) help ;;
        *) echo "Invalid flag ${flag}" && help ;;
    esac
done
shift $(( OPTIND-1 ))
[[ $# -eq 1 ]] && help
# main
compute_matrix_coverage "${@}"
