#!/bin/bash
set -euo pipefail
THREADS=$(nproc --all)
# MODE='cis'
MODE='gw'
COMMON_PARAMS="--mode ${MODE} --force"
# COMMON_PARAMS="--force --threads 12 --mode ${MODE}"
ICE_PARAMS="--in-memory --ignore-diags 2 --mad-max 5 --min-nnz 5 --min-count 0"
# SCALE_PARAMS="--max-percentile 10 --max-row-sum-err 0.05"
BALANCE_PARAM_SETS=(
    "ice   --name ${MODE}.ice   ${COMMON_PARAMS} ${ICE_PARAMS}"
    "vc    --name ${MODE}.vc    ${COMMON_PARAMS}"
    # "ice   --name ${MODE}.ice   ${COMMON_PARAMS} ${ICE_PARAMS}"
    # "scale --name ${MODE}.scale ${COMMON_PARAMS} ${SCALE_PARAMS}"
)
# GENOME_FASTA="/data/talkowski/Samples/16p_HiC/results/hg38.fasta.gz"
# GENOME_FASTA='/data/talkowski/tools/ref/Hi-c_ref/hg38.fasta'
# GENOME_FASTA='/home/sidreed/TalkowskiLab/Projects/HiC/remote.16p/qc.test/Hi-c_ref/hg38.fasta'
CHROMOSOMES=(
    "chr1"
    "chr2"
    "chr3"
    "chr4"
    "chr5"
    "chr6"
    "chr7"
    "chr8"
    "chr9"
    "chr10"
    "chr11"
    "chr12"
    "chr13"
    "chr14"
    "chr15"
    "chr16"
    "chr17"
    "chr18"
    "chr19"
    "chr20"
    "chr21"
    "chr22"
    "chrX"
)
balance_all_resolutions(){
    hic_samples=${@}
    # export HDF5_USE_FILE_LOCKING=FALSE  # necessary for hictk
    for param_set in "${BALANCE_PARAM_SETS[@]}"; do
    for sample_file in ${hic_samples[@]}; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        for uri in $(cooler ls ${sample_file}); do
            cooler balance --nproc 12 -f --name weight ${uri}
        done
        # normalizes all resolutions in the file
        echo "hictk balance $(echo "$param_set" | tr -s ' ') ${sample_file}"
        # hictk balance $(echo "$param_set" | tr -s ' ') ${sample_file}
    done
    done
}
dump_all_regions() {
    hic_samples=${@}
    for sample_file in ${hic_samples[@]}; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.mcool}"
        for uri in $(cooler ls ${sample_file}); do
        for region in ${REGIONS[@]}; do
            output_file="${OUTPUT_DIR}/${sample_ID}.${resolution}.${region}.txt"
            [[ -f "$output_file" ]] && continue
            basename "${output_file}"
            cooler dump \
                --nproc $THREADS \
                --range "${region}" \
                -t pixels \
                --join \
                --no-balance \
                --out "${output_file}" \
                "$uri"
        done
        done
    done
}
run_qc3c() {
    output_dir=$1
    hic_bams=${@:2}
    for sample_file in ${hic_bams[@]}; do 
        # samtools sort "${bam_file}" -o "${bam_file}"
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.lane1.hg38.0.bam}"
        echo $sample_ID
        qc3C bam \
            -t $(nproc --all) \
            -q 30 \
            --seed 9 \
            --max-obs 2000000 \
            --sample-rate 0.999 \
            -k arima \
            --fasta ${GENOME_FASTA} \
            --bam ${sample_file} \
            --output-path ${output_dir}/${sample_ID}
    done
}
digest_arima() {
    SIZES='/data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes'
    REFERENCE='/data/talkowski/tools/ref/Hi-c_ref/hg38.fasta'
    OUTPUT_DIR="/data/talkowski/Samples/16p_HiC/results/digestion"
    cooler digest $SIZES $REFERENCE DpnII > "${OUTPUT_DIR}/hg38.DpnII.restricted.bed"
    cooler digest $SIZES $REFERENCE HinfI > "${OUTPUT_DIR}/hg38.HinfI.restricted.bed"
    bedops --partition  "${OUTPUT_DIR}/hg38.DpnII.restricted.bed" "${OUTPUT_DIR}/hg38.HinfI.restricted.bed" >| "${OUTPUT_DIR}/hg36.ARIMA.restricted.bed"
}
digest() {
    SIZES='/data/talkowski/tools/ref/Hi-c_ref/hg38.reduced.chrom.sizes'
    REFERENCE='/data/talkowski/tools/ref/Hi-c_ref/hg38.fasta'
    DIGESTED_REFERENCE="/data/talkowski/Samples/16p_HiC/results/digestion/hg36.ARIMA.restricted.bed"
    SEED=9
    # FRAC=0.02 # ~ 1M reads getting sampled 
    output_dir=$1
    pairs_files=${@:2}
    for sample_file in ${pairs_files[@]}; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.nodups.pairs.gz}"
        output_file="${output_dir}/${sample_ID}.ARIMA.restricted.pairs"
        # normalizes all resolutions in the file
        pairtools restrict -f $DIGESTED_REFERENCE -o ${output_file} $sample_file
            # pairtools samples --seed $SEED $FRAC | \
    done
}
pair_stats() {
    output_dir=$1
    pairs_files=${@:2}
    for sample_file in ${pairs_files[@]}; do
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename "$sample_file")"
        sample_ID="${sample_ID%%.0.pairsam.gz}"
        sample_ID="${sample_ID/lane1./}"
        # normalizes all resolutions in the file
        # output_file="${output_dir}/${sample_ID}.pairs_stats.tsv"
        # pairtools stats -o ${output_file} $sample_file
        output_file="${output_dir}/${sample_ID}.scaling.tsv"
        pairtools scaling -o ${output_file} $sample_file
    done
}
# BALANCE mcool
mode="$1"
case $mode in 
    merge)
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
    balance) 
        balance_all_resolutions ${2} ${@:3}
        ;;
    dump) 
        dump_all_regions ${@:2}
        ;;
    'qc3')
        run_qc3c ${@:2}
        ;;
    'digest_chr')
        digest_arima 
        ;;
    digest)
        digest ${2} ${@:3}
        ;;
    stats)
        pair_stats ${@:2}
        ;;
    *)
        echo "Invalid mode: $mode" && exit 1
        ;;
esac
