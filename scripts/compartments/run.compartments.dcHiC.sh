#!/usr/bin/env bash
set -euo pipefail

###################################################
# Global Contstants
###################################################
GENOME_NAME="hg38"
DCHIC_DIR="$(pwd)/scripts/compartments"
PRE_PROCESS_SCRIPT="${DCHIC_DIR}/preprocess.dcHiC.py"
DCHIC_SCRIPT="${DCHIC_DIR}/run.dcHiC.R"
POST_PROCESS_SCRIPT="${DCHIC_DIR}/get_dcHiC_result.py"

###################################################
# Functions
###################################################
# Pre-process input for dcHiC
pre_process_samples() {
    output_dir="${1}"
    resolution="${2}"
    dcHiC_input_file="${3}"
    prefix1="${4}"
    cool_file1="${5}"
    prefix2="${6}"
    cool_file2="${7}"
    cd "${output_dir}"
    # pre-process first matrix
    matrix_file1="${output_dir}/${prefix1}.matrix"
    if ! [[ -e "${matrix_file1}" ]]; then
        python "${PRE_PROCESS_SCRIPT}"        \
            -genomeFile "${GENOME_SIZE_FILE}" \
            -res "${resolution}"              \
            -prefix "${prefix1}"              \
            -input cool                       \
            -file "${cool_file1}"
    fi
    # pre-process second matrix
    matrix_file2="${output_dir}/${prefix2}.matrix"
    if ! [[ -e "${matrix_file2}" ]]; then
        python "${PRE_PROCESS_SCRIPT}"        \
            -genomeFile "${GENOME_SIZE_FILE}" \
            -res "${resolution}"              \
            -prefix "${prefix2}"              \
            -input cool                       \
            -file "${cool_file2}"
    fi
# Generate input file for this comparison of conditions
    region_file1="${output_dir}/${prefix1}.abs.bed"
    region_file2="${output_dir}/${prefix2}.abs.bed"
cat > "${dcHiC_input_file}" <<EOF
${matrix_file1}	${region_file1}	${prefix1}	${prefix1}
${matrix_file2}	${region_file2}	${prefix2}	${prefix2}
EOF
    cd -
}
# Run dcHiC script from the original authors
run_dcHiC_pipeline() {
    output_dir="${1}"
    dcHiC_input_file="${2}"
    contact_type="${3}"
    pwd
    # Generate PCA loadings using cis contacts only
    Rscript "${DCHIC_SCRIPT}"        \
        --pcatype "${contact_type}"  \
        --dirovwt T                  \
        --cthread ${THREADS}         \
        --pthread ${THREADS}         \
        --file "${dcHiC_input_file}"
    # Pick which PCs to use for compartment assignment
    Rscript "${DCHIC_SCRIPT}"        \
        --pcatype select             \
        --genome "${GENOME_NAME}"    \
        --dirovwt T                  \
        --file "${dcHiC_input_file}"
    # Compute differential statistics between the 2 conditions
    Rscript "${DCHIC_DIR}"           \
        --pcatype analyze            \
        --diffdir "${output_dir}"    \
        --dirovwt T                  \
        --file "${dcHiC_input_file}"
    # Annotated subcompartments as well
    Rscript "${DCHIC_DIR}"           \
        --pcatype subcomp            \
        --diffdir "${output_dir}"    \
        --dirovwt T                  \
        --file "${dcHiC_input_file}"
}
# Clean up output data 
post_process_results() {
    output_dir="${1}"
    resolution="${2}"
    prefix1="${3}"
    prefix2="${4}"
    python "${POST_PROCESS_SCRIPT}"              \
        --genome_size_file "${GENOME_SIZE_FILE}" \
        --region_file "${REGION_FILE}"           \
        --res "${resolution}"                    \
        --prefix1 "${prefix1}"                   \
        --prefix2 "${prefix2}"                   \
        --datadir "${output_dir}/DifferentialResult/pcQnm" \
        --outfile "${output_dir}/${prefix1}-${prefix2}-dcHiC.results.tsv"
}
# run for all pairs of samples
main() {
    output_dir="${1}"
    resolution="${2}"
    contact_type="${3}"
    prefix1="${4}"
    cool_file1="${5}"
    prefix2="${6}"
    cool_file2="${7}"
    # set up args
    dcHiC_input_file="${output_dir}/${prefix1}-${prefix2}-dcHiC.input.txt"
    # region_file1="${output_dir}/${prefix1}-abs.bed"
    # region_file2="${output_dir}/${prefix2}-abs.bed"
    # run pipeline
    pre_process_samples "${output_dir}" "${resolution}" "${dcHiC_input_file}" "${prefix1}" "${cool_file1}" "${prefix2}" "${cool_file2}"
    run_dcHiC_pipeline "${output_dir}" "${dcHiC_input_file}" "${contact_type}"
    # post_process_results "${output_dir}" "${resolution}" "${prefix1}" "${prefix2}"
}

###################################################
# Handle args
###################################################
# Set up arguments
BASE_DIR="$(pwd)"
# DCHIC_DIR="${BASE_DIR}/scripts/compartments"
OUTPUT_DIR="${BASE_DIR}/results/compartments/example"
REF_DIR="${BASE_DIR}/reference.files/genome.reference"
GENOME_SIZE_FILE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
# REGION_FILE="${REF_DIR}/"
RESOLUTIONS=(100000 50000 25000 10000 5000)
RESOLUTION=${RESOLUTIONS[0]}
CONTACT_TYPES=('cis' 'trans')
CONTACT_TYPE=${CONTACT_TYPES[0]}
THREADS=8
# read cli args
while getopts :C:c:P:p:r:b:l:g:t:o:d:h opt; do
    case $opt in
        C) cool_file1=${OPTARG} ;;
        P) prefix1=${OPTARG} ;;
        c) cool_file2=${OPTARG} ;;
        p) prefix2=${OPTARG} ;;
        r) RESOLUTION=${OPTARG} ;;
        b) CONTACT_TYPE=${OPTARG} ;;
        # r) RESOLUTIONS=($(echo "${OPTARG}" | cut -d',' --output-delimiter=' ' -f1-)) ;;
        l) REGION_FILE=${OPTARG} ;;
        g) GENOME_SIZE_FILE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        o) OUTPUT_DIR=${OPTARG} ;;
        d) DCHIC_DIR=${OPTARG} ;;
        h) usage ;;
        *) usage ;;
    esac
done

output_dir="${OUTPUT_DIR}/contact.type_${CONTACT_TYPE}/resolution_${RESOLUTION}"
mkdir -p "${output_dir}"

main "${output_dir}" "${RESOLUTION}" "${CONTACT_TYPE}" "${prefix1}" "${cool_file1}" "${prefix2}" "${cool_file2}"
