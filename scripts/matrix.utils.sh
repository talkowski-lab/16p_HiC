#!/bin/bash
set -euo pipefail
CONDA_DIR="$HOME/miniforge3"
# Genomic reference files
REF_DIR="/data/talkowski/tools/ref/Hi_c_noalt"
GENOME_CHR_SIZES="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# Reference files for pairtools restrict using ARIMA HiC kit
REF_NAME=$(basename $GENOME_REFERENCE)
REF_NAME="${REF_DIR%%.fasta}"
DPNII_DIGESTION="${REF_DIR}/${REF_NAME}.DpnII.digested.bed"
HINFI_DIGESTION="${REF_DIR}/${REF_NAME}.HinfI.digested.bed"
ARIMA_DIGESTION="${REF_DIR}/hg38.ARIMA.restricted.bed"
SEED=9
THREADS="16" 
# Utils
activate_conda() {
    # activate conda env with specific tools for each task
    case "$1" in
        cooler)    env_name="cooltools" ;;
        cooltools) env_name="cooltools" ;;
        pairtools) env_name="dist2" ;;
        qc3C)      env_name="qc3c" ;;
        fanc)      env_name="fanc" ;;
        multiqc)   env_name="mqc" ;;
        *)      echo "Invalid conda env: $1" && exit 1 ;;
    esac
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}
dump_all_regions() {
    output_dir="$(readlink -e "${1}")"
    hic_samples=${@:2}
    mkdir -p "${output_dir}"
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
    output_dir="$1"
    output_dir="$(readlink -e "${1}")"
    region="$2"
    resolution="$3"
    hic_samples=${@:4}
    mkdir -p "${output_dir}"
    activate_conda 'fanc'
    for sample_file in ${hic_samples[@]}; do
        [[ ${sample_file} == *.mcool ]] && continue
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
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE DpnII >| "${DPNII_DIGESTION}"
    cooler digest $GENOME_CHR_SIZES $GENOME_REFERENCE HinfI >| "${HINFI_DIGESTION}"
    # "merge" the two digestions i.e. list all genome fragments has with a breakpoint at cut sites for 1 of any enzyme supplied
    bedops --partition "${DPNII_DIGESTION}" "${HINFI_DIGESTION}" >| "${ARIMA_DIGESTION}"
}
pairtools_restrict() {
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    mkdir -p "${output_dir}"
    activate_conda 'pairtools'
    for sample_file in ${pairs_files[@]}; do
        [[ ${sample_file} == *.nodups.pairs.gz ]] && continue
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
# QC Stats
make_multiqc_report() {
    report_file="${1}"
    input_dir=${2}
    multiqc \
        --no-ai \
        --no-data-dir \
        --clean-up \
        --force \
        --template default \
        --filename "${report_file}" \
        ${input_dir}
}
make_multiqc_reports() {
    distiller_output_dir="$(readlink -e ${1})"
    output_dir="${2}"
    report_dir="${output_dir}/multiqc.reports"
    mkdir -p "${report_dir}"
    activate_conda 'multiqc'
    make_multiqc_report \
        "${report_dir}/fastqc.multiqc" \
        ${distiller_output_dir}/fastqc/
    make_multiqc_report \
        "${report_dir}/fastp.multiqc" \
        ${distiller_output_dir}/mapped_parsed_sorted_chunks/
    make_multiqc_report \
        "${report_dir}/pairtools.multiqc" \
        ${distiller_output_dir}/pairs_library/
    make_multiqc_report \
        "${report_dir}/qc3C.multiqc" \
        ${distiller_output_dir}/qc3C/
    readlink -e ${report_dir}/*
# QC Stats
pairtools_stats() {
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    activate_conda 'pairtools'
    for sample_file in ${pairs_files[@]}; do
        [[ ${sample_file} == *.nodups.pairs.gz ]] && continue
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.hg38.nodups.pairs.gz}"
        pairtools scaling \
            --output "${output_dir}/${sample_ID}.scaling.tsv" \
            ${sample_file}
        pairtools stats \
            --bytile-dups \
            --with-chromsizes \
            --output "${output_dir}/${sample_ID}.stats.tsv" \
            ${sample_file}
    done
}
run_qc3c() {
    output_dir="$(readlink -e "${1}")"
    hic_bams=${@:2}
    activate_conda 'qc3C'
    mkdir -p "${output_dir}"
    echo ${output_dir}
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
        if [[ -e "${output_files_path}/report.qc3C.json" ]]; then
            echo "Skipping, cached results here: ${output_files_path}"
            continue
        fi
        # Run QC on 200,000 subsampled reads, all mapq > 30
        qc3C bam                        \
            --threads ${THREADS}        \
            --fasta ${GENOME_REFERENCE} \
            --library-kit arima         \
            --min-mapq 30               \
            --seed ${SEED}              \
            --max-obs 200000            \
            --sample-rate 0.5           \
            --bam "${sample_file}"      \
            --output-path "${output_files_path}"
    done
}
}
# Main, determine what to do and what input to expect
mode="$1"
case $mode in 
    "merge")
        MERGED_FILE="${2}"
        HIC_SAMPLES=${@:3}
        [[ -f "$MERGED_FILE" ]] && continue
        cooler merge "${MERGED_FILE}.cool" ${HIC_SAMPLES[@]}
        cooler zoomify \
            --out "${MERGED_FILE}.mcool" \
            --resolutions "${RESOLUTIONS}" \
            ${balance_flag} \
            "${MERGED_FILE}.cool"
        ;;
    merge_NIBPLWAPL) merge_NIPBLWAPL_matrices ;;
    merge_16p)       merge_16p_matrices ;;
    dump)            dump_all_regions ${@:2} ;;
    plot_triangle)   plot_fanc ${@:2} ;;
    digest_genome)   digest_genome_arima ;;
    restrict)        pairtools_restrict ${@:2} ;;
    qc3C)            run_qc3c ${@:2} ;;
    stats)           pairtools_stats ${@:2} ;;
    multiqcs)        make_multiqc_reports ${@:2} ;;
    *)               echo "Invalid mode: $mode" && exit 1 ;;
esac
