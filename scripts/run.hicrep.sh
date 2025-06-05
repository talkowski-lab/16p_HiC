#!/bin/bash
set -euo pipefail
# Fixed params
RESOLUTIONS=(500000 100000 50000 25000 10000)
get_smooting_param() {
    # get relevant smoothing param for each resolution
    # https://github.com/TaoYang-dev/hicrep
    resolution=$1
    case $resolution in 
        1000000) h=0  ;;
        500000)  h=1  ;;
        100000)  h=3  ;;
        50000)   h=4  ;;
        40000)   h=5  ;;
        25000)   h=10 ;;
        10000)   h=20 ;;
        5000)    h=40 ;;
    esac
    echo ${h}
}
MAX_WINDOW_SIZE=(100000 500000 1000000 5000000)
DOWNSAMPLE=(1 0)
# Unput files
OUTPUT_DIR="${1}"
mkdir -p "${OUTPUT_DIR}"
HIC_SAMPLES=${@:2}
# All combos of hicrep params
job_finished=0
jobs_skipped=0
for downsample in ${DOWNSAMPLE[@]}; do
for max_window_size in ${MAX_WINDOW_SIZE[@]}; do
for resolution in ${RESOLUTIONS[@]}; do
    h="$(get_smooting_param ${resolution})"
    # For each possible valid combination of parameters
    if [[ $downsample -eq 0 ]]; then
        downsample_arg="" 
        downsample_str="FALSE" 
    else
        downsample_arg="--bDownSample "
        downsample_str="TRUE"
    fi
    [[ $max_window_size -lt $resolution ]] && continue
    # For all combinations of samples
    for sample_file_i in ${HIC_SAMPLES[@]}; do
    for sample_file_j in ${HIC_SAMPLES[@]}; do
        # sample 1 info
        sample_file_i="$(readlink -e ${sample_file_i})"
        sample_ID_i="$(basename $sample_file_i)"
        sample_ID_i="${sample_ID_i%%.mcool}"
        # sample 2 info
        sample_file_j="$(readlink -e ${sample_file_j})"
        sample_ID_j="$(basename $sample_file_j)"
        sample_ID_j="${sample_ID_j%%.mcool}"
        # Check if results are redundant or already exist
        comparison_dir="resolution_${resolution}/h_${h}/is.downsampled_${downsample_str}/window.size_${max_window_size}"
        output_file="${OUTPUT_DIR}/${comparison_dir}/${sample_ID_i}-${sample_ID_j}-hicrep.txt"
        redundant_file="${OUTPUT_DIR}/${comparison_dir}/${sample_ID_j}-${sample_ID_i}-hicrep.txt"
        if [[ -e ${output_file} || -e ${redundant_file} ]]; then
            jobs_skipped=$(( $jobs_skipped+1 )) 
            echo ${jobs_skipped}
            continue 
        fi
        [[ ${sample_ID_i} == ${sample_ID_j} ]] && continue
        # Run hicrep
        mkdir -p "${OUTPUT_DIR}/${comparison_dir}"
        # echo $output_file
        cmd="hicrep ${downsample_arg}--dBPMax ${max_window_size} --binSize ${resolution} --h ${h} ${sample_file_i} ${sample_file_j} ${output_file}"
        echo ${cmd}
        ${cmd}
        # hicrep "${downsample_arg}"
        #     --h "${h}" 
        #     --dBPMax ${max_window_size} 
        #     --binSize ${resolution} 
        #     ${sample_file_i} 
        #     ${sample_file_j}
        #     ${output_file}
    done
    done
done
done
done
echo "Finished ${job_finished} HiCRep runs and skipped re-calcuating ${jobs_skipped} runs"
