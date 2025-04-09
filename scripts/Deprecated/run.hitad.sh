#!/bin/bash
# set -euo pipefail
# slurm params
# QUEUE="short"
THREADS=32
# MEM_GB=30
# hiTAD params
METHODS=(hiTAD)
# RESOLUTIONS=(500000 100000 50000 10000 5000)
RESOLUTIONS=(100000)
WEIGHTS=(weight)
WEIGHT_NAMES=(ICE)
# input files
OUTPUT_DIR="${1}"
LOG_DIR="${2}"
HIC_SAMPLES=${@:3}
# loop over all param combos
jobs_done=0
jobs_skipped=0
for method in ${METHODS[@]}; do
for resolution in ${RESOLUTIONS[@]}; do
for i in ${!WEIGHTS[@]}; do
    weight=${WEIGHTS[i]}
    weight_name=${WEIGHT_NAMES[i]}
    for sample_file in ${HIC_SAMPLES[@]}; do
        # sample info
        sample_file="$(readlink -e ${sample_file})"
        sample_ID="$(basename ${sample_file})"
        sample_ID="${sample_ID%%.mcool}"
        # log info
        param_dir="${method}/${weight_name}/${resolution}"
        job_name="$(echo $param_dir | sed -e 's/\//./g').${sample_ID}"
        log_file="${LOG_DIR}/${job_name}.log"
        # output files
        output_dir="${OUTPUT_DIR}/${param_dir}"
        output_TAD_file="${output_dir}/${sample_ID}-tads.txt"
        output_DI_file="${output_dir}/${sample_ID}-DIs.txt"
        # Dont submit job if results already exist
        [[ -a $output_TAD_file ]] && jobs_skipped=$(( jobs_skipped+1 )) || job_num=$(( job_num+1 ))
        [[ -a $output_TAD_file ]] && continue 
        mkdir -p "${output_dir}"
        echo $output_TAD_file
        # Run hitad
        domaincaller                      \
            --logFile ${log_file}         \
            --cpu-core $THREADS           \
            --weight-col ${weight}        \
            --output  ${output_TAD_file}  \
            --DI-output ${output_DI_file} \
            --uri "${sample_file}::resolutions/${resolution}"
    done
done
done
done
echo "Ran ${job_num} jobs, skipped ${jobs_skipped}"
