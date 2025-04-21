#!/bin/bash
set -euo pipefail
CONDA_DIR="$HOME/miniforge3"
# Genomic reference files
REF_DIR="/data/talkowski/tools/ref/Hi_c_noalt"
GENOME_CHR_SIZES="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
GENOME_REFERENCE="${REF_DIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
# Reference files for pairtools restrict using ARIMA HiC kit
REF_NAME=$(basename $GENOME_REFERENCE)
REF_NAME="${REF_NAME%%.fasta}"
DPNII_DIGESTION="${REF_DIR}/${REF_NAME}.DpnII.digested.bed"
HINFI_DIGESTION="${REF_DIR}/${REF_NAME}.HinfI.digested.bed"
SEED=9
THREADS="16" 
ARIMA_DIGESTION="${REF_DIR}/${REF_NAME}.ARIMA.digested.bed"
# HiC matrix relevant args
HIC_16p_RESULTS_DIR="/data/talkowski/Samples/16p_HiC"
HIC_NIPBLWAPL_RESULTS_DIR="/data/talkowski/Samples/WAPL_NIPBL/HiC"
MAPQ_FITERS=("mapq_30" "no_filter")
RESOLUTIONS="10000000,5000000,2500000,1000000,500000,250000,100000,50000,25000,10000,5000"
# Utils
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
    output_dir="$(readlink -e "${1}")"
    region="$2"
    resolution="$3"
    hic_samples=${@:4}
    mkdir -p "${output_dir}"
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
pairtools_restrict() {
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    mkdir -p "${output_dir}"
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
    # Merge contacts
    cool_file="${output_dir}/${merged_name}.hg38.${filter_status}.1000.cool"
    if ! [[ -f "${cool_file}" ]]; then
        cooler merge       \
            "${cool_file}" \
            ${hic_matrices[@]}
    fi
    # Bin + balance merged matrix at all specified resolutions
    mcool_file="${cool_file%%.cool}.mcool"
    if ! [[ -f "${mcool_file}" ]]; then
        cooler zoomify                     \
            --nproc ${THREADS}             \
            --resolutions "${RESOLUTIONS}" \
            --balance                      \
            --out "${mcool_file}"          \
            "${cool_file}"
    fi
}
merge_16p_matrices() {
    activate_conda 'cooler'
    cooler_dir="${HIC_16p_RESULTS_DIR}/results.NSC/coolers_library"
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
    cooler_dir="${HIC_NIPBLWAPL_RESULTS_DIR}/results.iN/coolers_library"
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
            ${cooler_dir}/WAP.DEL.{P1B5,P1E9,P2C4}.iN.HiC/*.${filter_name}.1000.cool
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
# QC Stats
pairtools_stats() {
    output_dir="$(readlink -e "${1}")"
    pairs_files=${@:2}
    mkdir -p "${output_dir}"
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
        # calculate pair stats
        stats_file="${output_dir}/${sample_ID}.stats.tsv"
        if ! [[ -f ${scale_file} ]]; then
            pairtools stats              \
                --bytile-dups            \
                --with-chromsizes        \
                --output "${stats_file}" \
                "${sample_file}"
        fi
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
    output_dir="$(readlink -e "${1}")"
    distiller_dir="$(readlink -e ${2})"
    mkdir -p "${output_dir}"
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
mode="$1"
case $mode in 
    merge_NIBPLWAPL) merge_NIPBLWAPL_matrices ;;
    merge_16p)       merge_16p_matrices ;;
    dump)            dump_all_regions ${@:2} ;;
    plot_triangle)   plot_fanc ${@:2} ;;
    digest_genome)   digest_genome_arima ;;
    restrict)        pairtools_restrict ${@:2} ;;
    qc3C)            run_qc3c ${@:2} ;;
    # stats)           pairtools_stats ${@:2} ;;
    multiqcs)        make_multiqc_reports ${@:2} ;;
    *)               echo "Invalid mode: $mode" && exit 1 ;;
esac
