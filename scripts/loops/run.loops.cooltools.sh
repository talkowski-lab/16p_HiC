#!/bin/bash
set -uo pipefail

###################################################
# Functions
###################################################
help() {
    echo "$(basename ${0}) [OPTIONS] sample{1..10}.mcool"
    exit 0
}

call_loops() {
    output_dir="$(readlink -e "${1}")"
    echo "Saving results in ${output_dir}"
    mkdir -p "${1}"
    for sample_file in "${@:2}"; do
        # Extract SampleID
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        sample_file="$(readlink -e ${sample_file})"
        # loop over all resolutions + normalizations + {cis,trans}
        # for uri in $(cooler ls "${sample_file}"); do 
            # resolution="$(echo "${uri}" | rev | cut -d '/' -f1 | rev)"
        for resolution in "${RESOLUTIONS[@]}"; do
            uri="${sample_file}::resolutions/${resolution}"
            for weight_name in "${!WEIGHTS[@]}"; do
                weight="${WEIGHTS[${weight_name}]}"
                if [[ ${weight} != '' ]]; then
                    weight_flag="--clr-weight-name ${weight}"
                else 
                    weight_flag=""
                fi
                for contact_type in "${CONTACT_TYPES[@]}"; do
                    if [[ ${contact_type} == 'cis' ]]; then
                        ct_args="--smooth --aggregate-smoothed"
                    else 
                        ct_args=""
                    fi
                    # name output file directory with params
                    param_dir="method_${METHOD}/type_${contact_type}/weight_${weight_name}/resolution_${resolution}"
                    # Caculate expected IF of each position
                    expected_contacts_file="${EXPECTED_CONTACTS_ROOT_DIR}/${param_dir}/${sample_ID}-expected.tsv"
                    if ! [[ -f "${expected_contacts_file}" ]]; then
                        echo "expected contacts  file doesnt exists, generating it"
                        mkdir -p "$(dirname "${expected_contacts_file}")"
                        cmd="cooltools expected-${contact_type} ${weight_flag} --nproc ${THREADS} ${ct_args} --output ${expected_contacts_file} ${uri}"
                        echo ${uri}
                        echo ${cmd}
                        ${cmd}
                    fi
                    output_file="${output_dir}/${param_dir}/${sample_ID}-dots.tsv"
                    # call dots with default args
                    [[ -f "${output_file}" ]] && continue
                    mkdir -p "$(dirname "${output_file}")"
                    cmd="cooltools dots ${weight_flag} --nproc ${THREADS} --output ${output_file} ${uri} ${expected_contacts_file}"
                    echo "${cmd}"
                    ${cmd}
            done
            done
        done
    done
}

###################################################
# Handle Arguments
###################################################
CONDA_DIR="${HOME}/miniforge3"
EXPECTED_CONTACTS_ROOT_DIR="./results/sample.QC/expected.coverage"
METHOD="cooltools"
# RESOLUTIONS=(100000 50000 25000 10000 5000)
RESOLUTIONS=(25000 10000 5000)
CONTACT_TYPES=('cis')
declare -rA WEIGHTS=([raw]='' [balanced]='weight')
THREADS=$(nproc)
# Handle CLI args
[[ $# -eq 0 ]] && echo "No Args" && exit 1
while getopts "a:t:m:h" flag; do
    case ${flag} in 
        m) METHOD="${OPTARG}" ;;
        e) EXPECTED_CONTACTS_ROOT_DIR="${OPTARG}" ;;
        a) CONDA_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) help ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))

###################################################
# Main
###################################################
# activate conda
source "${CONDA_DIR}/etc/profile.d/conda.sh"
conda activate "cooltools"
# call main fnc
call_loops "${@}"

