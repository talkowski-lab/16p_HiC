#!/usr/bin/env bash
set -euo pipefail

# Adapted from this script 
# https://github.com/pollardlab/contact_map_scoring/blob/d6f199759f72daff68b3672c66eef32c285df777/code/run_dcHiC.sh

# Functions

usage() {
    cat << EOF >&2
    Usage: run_dcHiC.sh [-C <cool_file1>] [-c <cool_file2>] [-r <region_file>] [-P <prefix1>]
                [-p <prefix2>] [-g <genome_size_file>] [-b <bin_size>] [-t <threads>] [-o <outdir>] [-d dchicdir]
    This script aims to run dcHiC from genome contact maps (.cool file) of two conditions.
    -C : The genome-wide contact map in .cool format from condition 1.
    -c : The genome-wide contact map in .cool format from condition 2.
    -r : The regions for map comparison.
    -P : The prefix of condition1.
    -p : The prefix of condition2.
    -g : The file of genome size.
    -b : bin size for contact map. Default is 2048.
    -t : The CPU threads for parallelized processing. Default is 1.
    -o : The output directory.
    -d : The directory of dcHiC.
    -h : Show usage help
EOF
    exit 1
}

run_dcHiC(){
    cool_file1="${1}"
    prefix1="${2}"
    cool_file2="${3}"
    prefix2="${4}"
    resolution="${5}"
    contact_type="${6}"
    output_dir="${7}"
    # Preprocess .cool files into matrices for dcHiC
    python "${PRE_PROCESS_SCRIPT}"        \
        -genomeFile "${GENOME_SIZE_FILE}" \
        -res "${resolution}"              \
        -prefix "${prefix1}"              \
        -input cool                       \
        -file "${cool_file1}"
    python "${PRE_PROCESS_SCRIPT}"        \
        -genomeFile "${GENOME_SIZE_FILE}" \
        -res "${resolution}"              \
        -prefix "${prefix2}"              \
        -input cool                       \
        -file "${cool_file2}"
    # Generate input file for this comparison of conditions
    input_file="${output_dir}/${prefix1}-${prefix2}-dcHiC.input.txt"
    region_file1="${output_dir}/${prefix1}-abs.bed"
    region_file2="${output_dir}/${prefix2}-abs.bed"
    cat > "${input_file}" <<EOF
${prefix1}.matrix	${region_file1}	${prefix1}	${prefix1}
${prefix2}.matrix	${region_file2}	${prefix2}	${prefix2}
EOF
    # Run dcHiC script from the original authors
    # Generate PCA loadings using cis contacts only
    Rscript "${DCHIC_SCRIPT}"       \
        --pcatype "${contact_type}" \
        --dirovwt T                 \
        --cthread ${THREADS}        \
        --pthread ${THREADS}        \
        --file "${input_file}"
    # Pick which PCs to use for compartment assignment
    Rscript "${DCHIC_SCRIPT}"     \
        --pcatype select          \
        --genome "${GENOME_NAME}" \
        --dirovwt T               \
        --file "${input_file}"
    # Compute differential statistics between the 2 conditions
    Rscript "${DCHIC_DIR}"        \
        --pcatype analyze         \
        --diffdir "${output_dir}" \
        --dirovwt T               \
        --file "${input_file}"
    # Annotated subcompartments as well
    Rscript "${DCHIC_DIR}"        \
        --pcatype subcomp         \
        --diffdir "${output_dir}" \
        --dirovwt T               \
        --file input.ES_NPC.txt 
    # Clean up output data 
    python "${POST_PROCESS_SCRIPT}"              \
        --genome_size_file "${GENOME_SIZE_FILE}" \
        --region_file "${REGION_FILE}"           \
        --res "${resolution}"                    \
        --prefix1 "${prefix1}"                   \
        --prefix2 "${prefix2}"                   \
        --datadir "${output_dir}/DifferentialResult/pcQnm" \
        --outfile "${output_dir}/${prefix1}-${prefix2}-dcHiC.results.tsv"
}

main() {
    cool_file1="${1}"
    prefix1="${2}"
    cool_file2="${3}"
    prefix2="${4}"
    for contact_type in "${CONTACT_TYPES[@]}"; do
    for resolution in "${RESOLUTIONS[@]}"; do
        param_dir="contact.type_${contact_type}/resolution_${resolution}"
        output_dir="${OUTPUT_DIR}/${param_dir}"
        run_dcHiC            \
           "${cool_file1}"   \
           "${prefix1}"      \
           "${cool_file2}"   \
           "${prefix2}"      \
           "${resolution}"   \
           "${contact_type}" \
           "${output_dir}"
    done
    done
}

# GLOBAL ARGS
GENOME_NAME="hg38"
# Default args
REF_DIR=""
REGION_FILE="${REF_DIR}/"
GENOME_SIZE_FILE="${REF_DIR}/"
CONTACT_TYPES=('cis')
# CONTACT_TYPES=('cis' 'trans')
RESOLUTIONS=""
THREADS=8
OUTPUT_DIR="./results/compartments/results/"
DCHIC_DIR="./scripts"
# Handle args
while getopts :C:c:r:P:p:g:b:t:o:d:h opt; do
    case $opt in
        C) cool_file1=${OPTARG};;
        c) cool_file2=${OPTARG};;
        r) region_file=${OPTARG};;
        P) prefix1=${OPTARG};;
        p) prefix2=${OPTARG};;
        g) GENOME_SIZE_FILE=${OPTARG};;
        b) RESOLUTIONS=${OPTARG};;
        t) THREADS=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        d) DCHIC_DIR=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done
# Set up arguments
PRE_PROCESS_SCRIPT="${DCHIC_DIR}/dcHiC.preprocess.py"
DCHIC_SCRIPT="${DCHIC_DIR}/dchicf.R"
# POST_PROCESS_SCRIPT="${DCHIC_DIR}/dcHiC.postprocess.R"
POST_PROCESS_SCRIPT="${DCHIC_DIR}/get_dcHiC_result.py"

[ -z "$prefix1" ]  && echo "Error! Please provide the prefix of condition1" && usage
[ -z "$prefix2" ]  && echo "Error! Please provide the prefix of condition2" && usage

if ! [[ -e "${cool_file1}" && -e "${cool_file2}" && -e "${region_file}" ]]; then
      echo "Error! Two genome-wide contact maps, a region file and a dcHiC input file are needed." && usage 
fi

cool_file1="$(readlink -e "${cool_file1}")"
cool_file2="$(readlink -e "${cool_file2}")"

main "${cool_file1}" "${cool_file2}"
