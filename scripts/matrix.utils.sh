#!/bin/bash
set -o pipefail

###################################################
# Fixed variables/lists 
###################################################
# 16p Information
GENOTYPES_16P=('WT' 'DEL' 'DUP')
CELLTYPES_16P=('iN' 'NSC')
# Edit information
PROJECT_EDITS=('NIPBL' 'WAPL' 'RAD21' 'CTCF')
# GENOTYPES_EDITS=('WT' 'DEL' 'BIALLELIC')
GENOTYPES_EDITS=('WT' 'DEL')
CELLTYPES_EDITS=('iN')
# Technical Args
SEED=9  # Random seed for qc3C
RESOLUTIONS=(100000 50000 25000 10000 5000)
# Genomic reference files
REF_DIR='/data/talkowski/tools/ref/Hi_c_noalt'
# GENOME_CHR_SIZES="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

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
    # hic_matrices=$(find "${cooler_dir}" -type f -name "${matrix_file_pattern}" | grep -vi 'merge' | paste -sd" ")
    hic_matrices=$(find "${cooler_dir}" -type f -name "${matrix_file_pattern}" | grep -vi 'merge' | grep -vE '16p.iN.WT.(FACS1|p44|p49).TR1' | paste -sd" ")
    echo ${hic_matrices}
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
            ${hic_matrices}
    fi
    # Bin + balance merged matrix at all specified resolutions
    merged_mcool_file="${merged_cool_file%%.cool}.mcool"
    if [[ -e ${merged_mcool_file} ]]; then
        echo "Skipping balancing, balanced merged file exists: ${merged_mcool_file}"
    else 
        echo "Balancing ${merged_cool_file} ..."
        resolutions="$(echo "${RESOLUTIONS[@]}" | sed -e 's/ /,/g')"
        cooler zoomify                     \
            --nproc "${THREADS}"           \
            --resolutions "${resolutions}" \
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
                # Merge all WTs across edits
                if [[ ${genotype} == 'WT' ]]; then
                    sample_group="${celltype}.${genotype}"
                    echo "All.${sample_group}"
                    echo '---------------'
                    merge_matrices        \
                        "${cooler_dir}"   \
                        "*.${sample_group}" \
                        "${read_filter}"  \
                        "All.${sample_group}.Merged.Merged"
                fi
            done
            # Merge CTCF.iN.BIALLELIC replicates explic
            sample_group="${edit}.${celltype}.BIALLELIC"
            matrix_file_pattern="${sample_group}.*.${read_filter}.1000.cool"
            merge_matrices        \
                "${cooler_dir}"   \
                "${sample_group}" \
                "${read_filter}"  \
                "${sample_group}.Merged.Merged"
            echo '===================================================================='
        done
    done
}

