#!/bin/bash
set -euo pipefail
# Technical Args
SEED=9  # Random seed for qc3C
# Genomic reference files
REF_DIR="/data/talkowski/tools/ref/Hi_c_noalt"
GENOME_CHR_SIZES="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# Reference files for pairtools restrict using ARIMA HiC kit
REF_NAME=$(basename $GENOME_REFERENCE)
REF_NAME="${REF_NAME%%.fasta}"
DPNII_DIGESTION="${REF_DIR}/${REF_NAME}.DpnII.digested.bed"
HINFI_DIGESTION="${REF_DIR}/${REF_NAME}.HinfI.digested.bed"
ARIMA_DIGESTION="${REF_DIR}/${REF_NAME}.ARIMA.digested.bed"
DDEI_DIGESTION="${REF_DIR}/${REF_NAME}.DdeI.digested.bed"
HIC3_DIGESTION="${REF_DIR}/${REF_NAME}.HIC3.digested.bed"
# MAPQ_FITERS=("mapq_30" "no_filter")
# CELLTYPES=("iN" "NSC")
CELLTYPES=("NSC")
MAPQ_FITERS=("mapq_30")
RESOLUTIONS="10000000,5000000,2500000,1000000,500000,250000,100000,50000,25000,10000,5000"
# 1000     #   1Kb
# 2000     #   2Kb
# 5000     #   5Kb
# 10000    #  10Kb
# 25000    #  25Kb
# 50000    #  50Kb
# 100000   # 100Kb
# 250000   # 250Kb
# 500000   # 500Kb
# 1000000  #   1Mb
# 2500000  # 2.5Mb
# 5000000  #   5Mb
# 10000000 #  10Mb
# Utils
help() {
    if [[ $# -eq 0 ]]; then
        echo "Usage: ${0} \${MODE}
modes: 
    coverage 
    digest_genome 
    restrict 
    merge_NIBPLWAPL 
    merge_16p 
    coverage 
    qc3C 
    multiqcs"
    else
        case ${1} in 
            dump)
                args="TODO"
                ;;
            plot_triangle)
                args="TODO"
                ;;
            digest_genome)
                args=""
                ;;
            restrict)
                args="\${OUTPUT_DIR} sample1.hg38.nodups.pairs.gz sample{2..N}.hg38.nodups.pairs.gz"
                ;;
            merge_NIBPLWAPL) 
                args="\${COOLER_DIR}"
                ;;
            merge_16p)
                args="\${COOLER_DIR}"
                ;;
            coverage)
                args="\${OUTPUT_DIR} sample1.mcool sample{2..N}.mcool"
                ;;
            qc3C)
                args="\${OUTPUT_DIR} \${ENZYME1} \${ENZYME2} sample1.lane1.hg38.0.bam sample{2..N}.lane1.hg38.0.bam"
                ;;
            multiqcs)
                args="\${OUTPUT_DIR} \${DISTILLER_OUTPUT_DIR}"
                ;;
            *) 
                echo "Invalid mode: $mode" && exit 1 
                ;;
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
    hic_samples=${@:2}
    activate_conda 'cooler'
    for sample_file in ${hic_samples[@]}; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        for uri in $(cooler ls ${sample_file}); do
        for chr in ${CHROMOSOMES[@]}; do
            output_file="${output_dir}/${sample_ID}.${resolution}.${region}.txt"
            [[ -f "$output_file" ]] && continue
            basename "${output_file}"
            cooler dump                \
                --nproc $THREADS       \
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
plot_fanc() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    region="$2"
    resolution="$3"
    hic_samples=${@:4}
    activate_conda 'fanc'
    for sample_file in ${hic_samples[@]}; do
        [[ ${sample_file} == *.mcool ]] || continue
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.*.mcool}"
        output_file="${output_dir}/${sample_ID}.${resolution}.${region}.heatmap.pdf"
        fancplot $region \
            -vv \
            --output ${output_file} \
            -p triangular \
            -u -l \
            --title "${sample_ID} Contacts on ${region%:*}" \
            ${sample_file}@${resolution}
        # fancplot chr16:0mb-96mb --output ./WT_VS_DEL.Merged.100K.heatmap.pdf -p split --title 'Merged WT vs Merged DEL Contacts chr16' ../coolers_library/16pWTNSCHIC_S5S6.hg38.mapq_30.1000.mcool::resolutions/100000 ../coolers_library/16pDELNSCHIC_S1S2.hg38.mapq_30.1000.mcool::resolutions/100000
done
}
# Restrict Fragment analysis
digest_genome_arima() {
    # The cooler^1 docs shows that to analyze a multi-enzyme digestion you can "partition" the two individual digestion, as bedops^2 does.
    ## 1: https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition
    ## 2: https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition
    activate_conda 'cooler'
    echo "Creating DpnII Ref at: ${DPNII_DIGESTION}"
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE DpnII >| "${DPNII_DIGESTION}"
    echo "Creating HinfI Ref at: ${HINFI_DIGESTION}"
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE HinfI >| "${HINFI_DIGESTION}"
    # "merge" the two digestions i.e. list all genome fragments has with a breakpoint at cut sites for 1 of any enzyme supplied
    echo "\"Merging\" digestions at: ${ARIMA_DIGESTION}"
    bedops --partition "${DPNII_DIGESTION}" "${HINFI_DIGESTION}" >| "${ARIMA_DIGESTION}"
}
digest_genome_hic3() {
    # The cooler^1 docs shows that to analyze a multi-enzyme digestion you can "partition" the two individual digestion, as bedops^2 does.
    ## 1: https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition
    ## 2: https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition
    activate_conda 'cooler'
    echo "Creating DpnII Ref at: ${DPNII_DIGESTION}"
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE DpnII >| "${DPNII_DIGESTION}"
    echo "Creating DdeI Ref at: ${DDEI_DIGESTION}"
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE HinfI >| "${DDEI_DIGESTION}"
    # "merge" the two digestions i.e. list all genome fragments has with a breakpoint at cut sites for 1 of any enzyme supplied
    echo "\"Merging\" digestions at: ${HIC3_DIGESTION}"
    bedops --partition "${DPNII_DIGESTION}" "${DDEI_DIGESTION}" >| "${HIC3_DIGESTION}"
}
pairtools_restrict() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    activate_conda 'pairtools'
    for sample_file in ${pairs_files[@]}; do
        [[ ${sample_file} == *.nodups.pairs.gz ]] || continue
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.hg38.nodups.pairs.gz}"
        output_file="${output_dir}/${sample_ID}.ARIMA.restricted.pairs"
        pairtools restrict               \
            --frags "${ARIMA_DIGESTION}" \
            --output ${output_file}      \
            ${sample_file}
    done
}
# Merging matrices 
merge_matrices() {
    # Arg list
    output_dir="${1}"
    merged_name="${2}"
    filter_status="${3}"
    hic_matrices=${@:4}
    output_dir="${output_dir}/${merged_name}"
    mkdir -p "${output_dir}"
    echo ${hic_matrices[@]}
    # Merge contacts
    cool_file="${output_dir}/${merged_name}.hg38.${filter_status}.1000.cool"
        cooler merge "${cool_file}" ${hic_matrices[@]}
    # Bin + balance merged matrix at all specified resolutions
    mcool_file="${cool_file%%.cool}.mcool"
        cooler zoomify                     \
            --nproc ${THREADS}             \
            --resolutions "${RESOLUTIONS}" \
            --balance                      \
            --out "${mcool_file}"          \
            "${cool_file}"
}
merge_16p_matrices() {
    cooler_dir="$(readlink -e "${1}")"
    activate_conda 'cooler'
    # Merge matrices with and without MAPQ filtering
    for filter_name in ${MAPQ_FITERS[@]}; do 
        # WTs
        merge_matrices            \
            "${cooler_dir}"       \
            16p.WT.Merged.NSC.HiC \
            ${filter_name}        \
            ${cooler_dir}/16p.WT.{p46,FACS1,p49}.NSC.HiC/*.${filter_name}.1000.cool
        # DELs
        merge_matrices             \
            "${cooler_dir}"        \
            16p.DEL.Merged.NSC.HiC \
            ${filter_name}         \
            ${cooler_dir}/16p.DEL.{A3,B8,H10}.NSC.HiC/*.${filter_name}.1000.cool
        # DUPs
        merge_matrices             \
            "${cooler_dir}"        \
            16p.DUP.Merged.NSC.HiC \
            ${filter_name}         \
            ${cooler_dir}/16p.DUP.{C5,D12,G7}.NSC.HiC/*.${filter_name}.1000.cool
    done
}
merge_NIPBLWAPL_matrices() {
    activate_conda 'cooler'
    cooler_dir="$(readlink -e "${1}")"
    # Merge matrices with and without MAPQ filtering
    for filter_name in ${MAPQ_FITERS[@]}; do 
        # NIBPL WTs 
        merge_matrices             \
            "${cooler_dir}"        \
            NIBPL.WT.Merged.iN.HiC \
            ${filter_name}         \
            ${cooler_dir}/NIPBL.WT.{P1A1,P2B6,P2F4}.iN.HiC/*.${filter_name}.1000.cool
        # NIBPL DELs 
        merge_matrices              \
            "${cooler_dir}"         \
            NIBPL.DEL.Merged.iN.HiC \
            ${filter_name}          \
            ${cooler_dir}/NIPBL.DEL.{P1H10,P2E1,P2E6}.iN.HiC/*.${filter_name}.1000.cool
        # WAPL WTs
        merge_matrices            \
            "${cooler_dir}"       \
            WAPL.WT.Merged.iN.HiC \
            ${filter_name}        \
            ${cooler_dir}/WAPL.WT.{P1B10,P1E5,P2F2}.iN.HiC/*.${filter_name}.1000.cool
        # WAPL DELs
        merge_matrices             \
            "${cooler_dir}"        \
            WAPL.DEL.Merged.iN.HiC \
            ${filter_name}         \
            ${cooler_dir}/WAPL.DEL.{P1B5,P1E9,P2C4}.iN.HiC/*.${filter_name}.1000.cool
        # All WTs 
        merge_matrices           \
            "${cooler_dir}"      \
            All.WT.Merged.iN.HiC \
            ${filter_name}       \
            ${cooler_dir}/NIPBL.WT.{P1A1,P2B6,P2F4}.iN.HiC/*.${filter_name}.1000.cool ${cooler_dir}/WAPL.WT.{P1B10,P1E5,P2F2}.iN.HiC/*.${filter_name}.1000.cool
        # # All DELs 
        # merge_matrices             \
        #     "${cooler_dir}"        \
        #     \
        #     ${filter_name}         \
        #     ${cooler_dir}//*.${filter_name}.1000.cool
    done
}
# Merging stats files
merge_stats() {
    # Arg list
    output_dir="${1}"
    merged_name="${2}"
    stats_files="${3}"
    output_dir="${output_dir}/${merged_name}"
    mkdir -p "${output_dir}"
    ls ${stats_files}
    # Merge contacts
    # output_pairs_file="${output_dir}/${merged_name}.hg38.nodups.pairs.gz"
    # pairtools merge \
    #     --nproc ${THREADS} \
    #     --output ${output_pairs_file} \
    #     ${pairs_files[@]}
    output_stats_file="${output_dir}/${merged_name}.hg38.dedup.stats"
    pairtools stats \
        --merge \
        --output ${output_stats_file} \
        ${stats_files}
        
}
merge_16p_stats() {
    pairs_dir="$(readlink -e "${1}")"
    activate_conda 'cooler'
    # Merge pairs with and without MAPQ filtering
    for celltype in ${CELLTYPES[@]}; do 
        # WTs
        merge_stats \
            "${pairs_dir}" \
            16p.WT.Merged.${celltype}.HiC \
            "${pairs_dir}/**/16p.WT.*.${celltype}.HiC.hg38.dedup.stats"
            # $(find ${pairs_dir} -type f -name "16p.WT.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # DELs
        merge_stats \
            "${pairs_dir}" \
            "16p.DEL.Merged.${celltype}.HiC" \
            "${pairs_dir}/**/16p.DEL.*.${celltype}.HiC.hg38.dedup.stats"
            # $(find ${pairs_dir} -type f -name "16p.DEL.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # DUPs
        merge_stats \
            "${pairs_dir}" \
            "16p.DUP.Merged.${celltype}.HiC" \
            "${pairs_dir}/**/16p.DUP.*.${celltype}.HiC.hg38.dedup.stats"
            # $(find ${pairs_dir} -type f -name "16p.DUP.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
    done
}
merge_NIPBLWAPL_stats() {
    activate_conda 'cooler'
    pairs_dir="$(readlink -e "${1}")"
    # Merge pairs with and without MAPQ filtering
    for celltype in ${CELLTYPES[@]}; do 
        # NIBPL WTs 
        merge_pairs \
            "${pairs_dir}" \
            NIBPL.WT.Merged.${celltype}.HiC \
            $(find ${pairs_dir} -type f -name "NIPBL.WT.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # NIBPL DELs 
        merge_pairs \
            "${pairs_dir}" \
            NIBPL.DEL.Merged.${celltype}.HiC \
            $(find ${pairs_dir} -type f -name "NIPBL.DEL.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # WAPL WTs
        merge_pairs \
            "${pairs_dir}" \
            WAPL.WT.Merged.${celltype}.HiC \
            $(find ${pairs_dir} -type f -name "WAPL.WT.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # WAPL DELs
        merge_pairs             \
            "${pairs_dir}"        \
            WAPL.DEL.Merged.${celltype}.HiC \
            $(find ${pairs_dir} -type f -name "WAPL.DEL.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # All WTs 
        merge_pairs \
            "${pairs_dir}" \
            All.WT.Merged.${celltype}.HiC \
            $(find ${pairs_dir} -type f -name "*.WT.*.${celltype}.HiC.hg38.dedup.stats" | grep -v ".Merged.")
        # # All DELs 
        # merge_pairs             \
        #     "${pairs_dir}"        \
        #     \
        #     ${filter_name}         \
        #     ${pairs_dir}//*.${filter_name}.1000.cool
    done
}
# QC Stats
matrix_coverage() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    hic_matrices=${@:2}
    activate_conda 'cooltools'
    for sample_file in ${hic_matrices[@]}; do
        [[ ${sample_file} == *.mcool ]] || continue
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        sample_file="$(readlink -e "${sample_file}")"
        for uri in $(cooler ls ${sample_file}); do 
            resolution="$(echo "${uri}" | rev | cut -d '/' -f1 | rev)"
            raw_output_dir="${output_dir}/weight_raw/resolution_${resolution}"
            mkdir -p ${raw_output_dir}
            raw_output_file="${raw_output_dir}/${sample_ID}-coverage.tsv"
            echo ${uri}
            if ! [[ -e ${raw_output_file} ]]; then
                echo ${raw_output_file}
                cooltools coverage \
                    --nproc ${THREADS} \
                    --output ${raw_output_file} \
                    ${uri}
            fi
            balanced_output_dir="${output_dir}/weight_balanced/resolution_${resolution}"
            mkdir -p ${balanced_output_dir}
            balanced_output_file="${balanced_output_dir}/${sample_ID}-coverage.tsv"
            if ! [[ -e ${balanced_output_file} ]]; then
                echo ${balanced_output_file}
                cooltools coverage \
                    --nproc ${THREADS} \
                    --clr_weight_name 'weight' \
                    --output ${balanced_output_file} \
                    ${uri}
            fi
        done
    done
}
pairtools_stats() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    activate_conda 'pairtools'
    for sample_file in ${pairs_files[@]}; do
        [[ ${sample_file} == *.nodups.pairs.gz ]] || continue
        sample_file="$(readlink -e "${sample_file}")"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.hg38.nodups.pairs.gz}"
        echo ${sample_ID}
        # calculate dist/freq scaling
        scale_file="${output_dir}/${sample_ID}.scaling.tsv"
        if ! [[ -f ${scale_file} ]]; then
            pairtools scaling            \
                --output "${scale_file}" \
                "${sample_file}"
        fi
        # # calculate pair stats
        # stats_file="${output_dir}/${sample_ID}.stats.tsv"
        # if ! [[ -f ${scale_file} ]]; then
        #     pairtools stats              \
        #         --bytile-dups            \
        #         --with-chromsizes        \
        #         --output "${stats_file}" \
        #         "${sample_file}"
        # fi
    done
}
run_qc3c() {
    mkdir -p "${1}"
    output_dir="$(readlink -e "${1}")"
    enzyme1=${2}
    enzyme2=${3}
    hic_bams=${@:4}
    activate_conda 'qc3C'
    for sample_file in ${hic_bams[@]}; do 
        echo $sample_file
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.lane1.hg38.0.bam}"
        echo $sample_ID
        # Sort bamfile if not already sorted
        if ! $(samtools view -H "${sample_file}" | grep -q "SO:queryname"); then
            echo "Sorting bam file..."
            samtools sort            \
                --threads ${THREADS} \
                -n                   \
                -o "${sample_file}"  \
                "${sample_file}"
        fi
        # Skip if qc3C result already exists
        output_files_path="${output_dir}/${sample_ID}"
        echo $output_files_path
        # if [[ -e "${output_files_path}/report.qc3C.json" ]]; then
        #     echo "Skipping, cached results here: ${output_files_path}"
        #     continue
        # fi
        # Run QC on 200,000 subsampled reads, all mapq > 30
        qc3C bam                        \
            --threads ${THREADS}        \
            --fasta ${GENOME_REFERENCE} \
            --enzyme ${enzyme1}         \
            --enzyme ${enzyme2}         \
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
    echo ${tmp_file}
    echo find ${input_dir} -maxdepth 99 -type f -name "*${file_ext}"
    find ${input_dir} -maxdepth 99 -type f -name "*${file_ext}" >| ${tmp_file}
    multiqc                         \
        --no-ai                     \
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
    distiller_dir="$(readlink -e ${2})"
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
        # "${distiller_dir}/sample.QC/pairtools.stats/" \
    make_multiqc_report                               \
        "${output_dir}/pairtools.multiqc"             \
        "${distiller_dir}/pairs_library/" \
        ".dedup.stats"
    make_multiqc_report                    \
        "${output_dir}/qc3C.multiqc"       \
        "${distiller_dir}/sample.QC/qc3C/" \
        "report.qc3C.json"
    # list all generated reports
    readlink -e ${output_dir}/*
}
# Main, determine what to do and what input to expect
main() {
    mode="$1"
    case $mode in 
        merge_NIPBLWAPL_stats) 
            merge_NIPBLWAPL_stats ${2} 
            ;;
        merge_NIPBLWAPL_matrices) 
            merge_NIPBLWAPL_matrices ${2} 
            ;;
        merge_16p_stats) 
            merge_16p_stats ${2} 
            ;;
        merge_16p_matrices) 
            merge_16p_matrices ${2}
            ;;
        dump)
            dump_all_regions ${@:2}
            ;;
        plot_triangle)
            plot_fanc ${@:2} 
            ;;
        digest_genome)   
            digest_genome_arima
            digest_genome_hic3
            ;;
        restrict) 
            pairtools_restrict ${@:2} 
            ;;
        coverage) 
            matrix_coverage ${@:2} 
            ;;
        qc3C) 
            run_qc3c ${@:2} 
            ;;
        # stats)
        #     pairtools_stats ${@:2} 
        #     ;;
        multiqcs)
            make_multiqc_reports ${@:2} 
            ;;
        *) 
            echo "Invalid mode: $mode" && exit 1 
            ;;
    esac
}
# Default script arguments
CONDA_DIR="$HOME/miniforge3"
THREADS="32" 
# Handle CLI args
[[ $? -ne 0 ]] && echo "No Args" && exit 1
VALID_ARGS=$(getopt -o ha:t: --long help,anaconda-dir,threads -- "$@")
eval set -- "$VALID_ARGS"
while [ : ]; do
    case "$1" in
        -a|--anaconda-dir)
            CONDA_DIR="${2}"
            shift 2
            ;;
        -t|--threads)
            THREADS="${2}"
            shift 2
            ;;
        -h|--help) 
            shift 1 && help
            ;;
        --)
            shift 
            break
            ;;
    esac
done
if [[ $# -eq 1 ]]; then
    help ${1}
else 
    main ${@}
fi
