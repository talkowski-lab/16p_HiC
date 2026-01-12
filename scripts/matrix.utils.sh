#!/bin/bash
# set -euo pipefail
set -uo pipefail

###################################################
# Fixed variables/lists 
###################################################
# 16p Information
GENOTYPES_16P=('WT' 'DEL' 'DUP')
CELLTYPES_16P=('iN' 'NSC')
# CLONEIDS_16P=('A12' 'A3' 'B8' 'C5' 'D12' 'D9' 'FACS1' 'G7' 'H10' 'p44' 'p46' 'p49')
# Edit information
PROJECT_EDITS=('NIPBL' 'WAPL' 'RAD21')
GENOTYPES_EDITS=('WT' 'DEL')
CELLTYPES_EDITS=('iN')
# Technical Args
SEED=9  # Random seed for qc3C
MAPQ_FITERS=('mapq_30')
# RESOLUTIONS='10000000,5000000,2500000,1000000,500000,250000,100000,50000,25000,10000,5000'
# RESOLUTIONS=(10000000 5000000 2500000 1000000 500000 250000 100000 50000 25000 10000 5000)
RESOLUTIONS=(100000 50000 25000 10000 5000)
CONTACT_TYPES=('cis' 'trans')
WEIGHT_NAMES=('raw' 'balanced')
declare -rA WEIGHTS=([raw]='' [balanced]='weight')
# Genomic reference files
REF_DIR='/data/talkowski/tools/ref/Hi_c_noalt'
# GENOME_CHR_SIZES="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

###################################################
# Utils
###################################################
help() {
    if [[ -v $1 ]]; then
echo "Usage: ${0} \${MODE}
modes: 
    merge_Cohesin 
    merge_16p 
    qc3C 
    multiqcs
    coverage 
    expected"
        exit 0
    else 
        case ${1} in 
            merge_Cohesin) args="\${COOLER_DIR}" ;;
            merge_16p)     args="\${COOLER_DIR}" ;;
            qc3C)          args="\${OUTPUT_DIR} \${ENZYME1} \${ENZYME2} sample{1..N}.lane1.hg38.0.bam" ;;
            multiqcs)      args="\${OUTPUT_DIR} \${DISTILLER_OUTPUT_DIR}" ;;
            coverage)      args="\${OUTPUT_DIR} sample{1..N}.mcool" ;;
            expected)      args="\${OUTPUT_DIR} sample{1..N}.mcool" ;;
            *) echo "Invalid mode: $mode" && exit 1 ;;
        esac
        echo "Usage: ${0} ${1} ${args}"
    fi
    exit 0
}

activate_conda() {
    # activate conda env with specific tools for each task
    case "$1" in
        cooler)    env_name="cooltools" ;;
        cooltools) env_name="cooltools" ;;
        pairtools) env_name="pairtools" ;;
        qc3C)      env_name="qc3c" ;;
        fanc)      env_name="fanc" ;;
        multiqc)   env_name="mqc" ;;
        *)      echo "Invalid conda env: $1" && exit 1 ;;
    esac
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

dump_all_regions() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    activate_conda 'cooler'
    # hic_samples=${@:2}
    # for sample_file in ${hic_samples[@]}; do
    for sample_file in "${@:2}"; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        for uri in $(cooler ls "${sample_file}"); do
        for chr in "${CHROMOSOMES[@]}"; do
            param_dir="${output_dir}/weight_${weight_name}/resolution_${resolution}/chr_${chr}"
            mkdir -p "${param_dir}"
            output_file="${param_dir}/${sample_ID}.matrix.txt.gz"
            [[ -f "$output_file" ]] && continue
            basename "${output_file}"
            cooler dump                \
                --nproc "${THREADS}"   \
                --range "${chr}"       \
                --table pixels         \
                --join                 \
                --no-balance           \
                --out "${output_file}" \
                "$uri"
        done
        done
    done
}

###################################################
# QC Stats
###################################################
run_qc3c() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    enzyme1=${2}
    enzyme2=${3}
    activate_conda 'qc3C'
    for sample_file in "${@:4}"; do 
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.lane1.hg38.0.bam}"
        # echo "${sample_ID}" && continue
        # Sort bamfile if not already sorted
        if ! samtools view -H "${sample_file}" | grep -q "SO:queryname"; then
            echo "Sorting bam file..."
            samtools sort              \
                --threads "${THREADS}" \
                -n                     \
                -o "${sample_file}"    \
                "${sample_file}"
        fi
        # Skip if qc3C result already exists
        output_files_path="${output_dir}/${sample_ID}"
        if [[ -f "${output_files_path}/report.qc3C.json" ]]; then
            echo "Skipping ${sample_ID} cuz results exist"
            continue
        fi
        # Run QC on 200,000 subsampled reads, all mapq > 30
        qc3C bam                        \
            --threads "${THREADS}"      \
            --fasta ${GENOME_REFERENCE} \
            --enzyme "${enzyme1}"       \
            --enzyme "${enzyme2}"       \
            --min-mapq 30               \
            --seed ${SEED}              \
            --max-obs 200000            \
            --sample-rate 0.5           \
            --bam "${sample_file}"      \
            --output-path "${output_files_path}"
    done
}

make_multiqc_report() {
    report_file="${1}"
    input_dir="${2}"
    file_ext="${3}"
    tmp_file=$(mktemp)
    # echo ${tmp_file}
    # echo find ${input_dir} -maxdepth 99 -type f -name "*${file_ext}"
    find "${input_dir}" -maxdepth 99 -type f -name "*${file_ext}" >| "${tmp_file}"
        # --no-ai                     \
    multiqc                         \
        --no-data-dir               \
        --clean-up                  \
        --force                     \
        --template default          \
        --filename "${report_file}" \
        --file-list "${tmp_file}"
}

make_multiqc_reports() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    distiller_dir="$(readlink -e "${2}")"
    activate_conda 'multiqc'
    # make each multi-sample report by datatype
    make_multiqc_report                \
        "${output_dir}/fastqc.multiqc" \
        "${distiller_dir}/fastqc/"     \
        "_fastqc.zip"
    make_multiqc_report                                 \
        "${output_dir}/fastp.multiqc"                   \
        "${distiller_dir}/mapped_parsed_sorted_chunks/" \
        ".fastp.json"
    make_multiqc_report                   \
        "${output_dir}/pairtools.multiqc" \
        "${distiller_dir}/pairs_library/" \
        ".dedup.stats"
    make_multiqc_report                    \
        "${output_dir}/qc3C.multiqc"       \
        "${distiller_dir}/sample.QC/qc3C/" \
        "report.qc3C.json"
    # list all generated reports
    find "${output_dir}" -type f -exec readlink -e {} \;
}

###################################################
# Merging matrices 
###################################################
merge_matrices() {
    # Arg list
    output_dir="$(readlink -e "${1}")"
    mkdir -p "${output_dir}"
    sample_group="${2}"
    read_filter="${3}"
    merged_sample_name="${4}"
    # list all samples to merge, ignore existing merged matrices
    matrix_file_pattern="${sample_group}.*.${read_filter}.1000.cool"
    hic_matrices=$(find "${cooler_dir}" -type f -name "${matrix_file_pattern}" | grep -vi 'merge' | paste -sd" ")
    # hic_matrices=$(find "${cooler_dir}" -type f -name "${matrix_file_pattern}" | grep -vi 'merge' | grep -vE '16p.iN.WT.(FACS1|p44|p49).TR1' | paste -sd" ")
    # Merge contacts
    sample_group_dir="${output_dir}/${merged_sample_name}"
    mkdir -p "${sample_group_dir}"
    merged_cool_file="${sample_group_dir}/${merged_sample_name}.hg38.${read_filter}.1000.cool"
    if [[ -e ${merged_cool_file} ]]; then
        echo 'Skipping merge, merged file exists'
    else 
        echo "Merging ${sample_group} samples into ${merged_cool_file}"
        cooler merge              \
            "${merged_cool_file}" \
            "${hic_matrices[@]}"
    fi
    # Bin + balance merged matrix at all specified resolutions
    merged_mcool_file="${merged_cool_file%%.cool}.mcool"
    if [[ -e ${merged_mcool_file} ]]; then
        echo "Skipping balancing, balanced merged file exists: ${merged_mcool_file}"
    else 
        echo "Balancing ${merged_cool_file} ..."
        cooler zoomify                     \
            --nproc "${THREADS}"           \
            --resolutions "$(echo "${RESOLUTIONS[@]}" | paste -sd',')" \
            --balance                      \
            --out "${merged_mcool_file}"   \
            "${merged_cool_file}"
    fi
}

merge_16p_matrices() {
    cooler_dir="$(readlink -e "${1}")"
    activate_conda 'cooler'
    # echo ${cooler_dir}
    # Merge matrices with and without MAPQ filtering
    for read_filter in "${MAPQ_FITERS[@]}"; do 
    for celltype in "${CELLTYPES_16P[@]}"; do 
    for genotype in "${GENOTYPES_16P[@]}"; do 
        # # merge across all technical replicates
        # for cloneID in ${CLONEIDS_16P[@]}; do 
        #     sample_group="16p.${celltype}.${genotype}.${cloneID}"
        #     echo "${sample_group}"
        #     echo '---------------'
        #     merge_matrices        \
        #         "${cooler_dir}"   \
        #         "${sample_group}" \
        #         "${read_filter}"  \
        #         "${sample_group}.Merged"
        # done
        # merge across all biological + technical replicates
        sample_group="16p.${celltype}.${genotype}"
        echo "${sample_group}"
        echo '--------------------------------------------------------------------'
        merge_matrices        \
            "${cooler_dir}"   \
            "${sample_group}" \
            "${read_filter}"  \
            "${sample_group}.Merged.Merged"
        echo '===================================================================='
    done
    done
    done
}

merge_Cohesin_matrices() {
    activate_conda 'cooler'
    cooler_dir="$(readlink -e "${1}")"
    # For all groups of matrices
    for read_filter in "${MAPQ_FITERS[@]}"; do 
    for celltype in "${CELLTYPES_EDITS[@]}"; do 
    for genotype in "${GENOTYPES_EDITS[@]}"; do 
        # merge across all biological + technical replicates PER Edit
        for edit in "${PROJECT_EDITS[@]}"; do 
            sample_group="${edit}.${celltype}.${genotype}"
            matrix_file_pattern="${sample_group}.*.${read_filter}.1000.cool"
            echo "${sample_group}"
            echo '---------------------------------------------'
            merge_matrices        \
                "${cooler_dir}"   \
                "${sample_group}" \
                "${read_filter}"  \
                "${sample_group}.Merged.Merged"
            echo '============================================='
        done
        # Merge across edits per genotype & celltype ACROSS Edits
        sample_group="${celltype}.${genotype}"
        echo "${sample_group}"
        echo '---------------'
        merge_matrices        \
            "${cooler_dir}"   \
            "*.${sample_group}" \
            "${read_filter}"  \
            "All.${sample_group}.Merged.Merged"
    done
    done
    done
}

###################################################
# Compute stuff
###################################################
matrix_balance() {
    activate_conda 'cooltools'
    for sample_file in "${@}"; do
        # Extract SampleID
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        sample_file="$(readlink -e "${sample_file}")"
        # for uri in $(cooler ls "${sample_file}"); do 
        for resolution in "${RESOLUTIONS[@]}"; do
            uri="${sample_file}::resolutions/${resolution}"
            echo "${uri}"
            # Caculate balancing weights
            cooler balance           \
                --nproc "${THREADS}" \
                "${uri}"
            done
    done
}

matrix_coverage() {
    output_dir="$(readlink -e "${1}")"
    mkdir -p "${1}"
    echo "Saving results in ${output_dir}"
    activate_conda 'cooltools'
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
            for weight_name in "${WEIGHT_NAMES[@]}"; do
                weight="${WEIGHTS[${weight_name}]}"
                if [[ ${weight} != '' ]]; then
                    weight_flag="--clr_weight_name ${weight} "
                else 
                    weight_flag=""
                fi
                # name output file directory with params
                param_dir="${output_dir}/weight_${weight_name}/resolution_${resolution}"
                mkdir -p "${param_dir}"
                output_file="${param_dir}/${sample_ID}-coverage.tsv"
                # echo "${output_file}"
                [[ -f "${output_file}" ]] && continue # skip if output file exists
                # Caculate total IF of each bin
                echo "${output_file}"
                cooltools coverage ${weight_flag}--nproc "${THREADS}" \
                    --output "${output_file}" \
                    "${uri}"
            done
        done
    done
}

matrix_expected() {
    output_dir="$(readlink -e "${1}")"
    echo "Saving results in ${output_dir}"
    mkdir -p "${1}"
    activate_conda 'cooler'
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
            for weight_name in "${WEIGHT_NAMES[@]}"; do
                weight="${WEIGHTS[${weight_name}]}"
                if [[ ${weight} != '' ]]; then
                    weight_flag="--clr-weight-name ${weight} "
                else 
                    weight_flag=""
                fi
                for contact_type in "${CONTACT_TYPES[@]}"; do
                    if [[ ${contact_type} == 'cis' ]]; then
                        ct_args="--smooth --aggregate-smoothed "
                    else 
                        ct_args=""
                    fi
                    # name output file directory with params
                    param_dir="${output_dir}/type_${contact_type}/weight_${weight_name}/resolution_${resolution}"
                    mkdir -p "${param_dir}"
                    output_file="${param_dir}/${sample_ID}-expected.tsv"
                    [[ -f "${output_file}" ]] && continue
                    echo "${output_file}"
                    # Caculate expected IF of each position
                    cmd="cooltools expected-${contact_type} ${weight_flag}--nproc ${THREADS} ${ct_args} --output ${output_file} ${uri}"
                    ${cmd}
                    # cooltools "expected-${contact_type}" ${weight_flag}--nproc "${THREADS}" \
                    #     --smooth \
                    #     --aggregate-smoothed \
                    #     --output "${output_file}" \
                    #     "${uri}"
            done
            done
        done
    done
}

###################################################
# Main
###################################################
main() {
    mode="$1"
    echo "${mode}"
    case ${mode} in 
        merge_Cohesin) merge_Cohesin_matrices "${2}" ;;
        merge_16p)     merge_16p_matrices "${2}" ;;
        qc3C)          run_qc3c "${@:2}" ;;
        multiqcs)      make_multiqc_reports "${@:2}" ;;
        dump)          dump_all_regions "${@:2}" ;;
        balance)       matrix_balance "${@:2}" ;;
        coverage)      matrix_coverage "${@:2}" ;;
        expected)      matrix_expected "${@:2}" ;;
        *) echo "Invalid mode: ${mode}" && exit 1 ;;
    esac
}

###################################################
# Handle Arguments
###################################################
CONDA_DIR="${HOME}/miniforge3"
THREADS="8" 
# Handle CLI args
[[ $# -eq 0 ]] && echo "No Args" && exit 1
while getopts "a:t:h" flag; do
    case ${flag} in 
        a) CONDA_DIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) help ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))
[[ $# -eq 1 ]] && help "${1}"
main "${@}"
