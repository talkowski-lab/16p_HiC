#!/bin/bash
set -euo pipefail

RESOLUTIONS=(500000 100000 50000 10000)
HS=(1 3 5 20)
MAX_WINDOW_SIZE=(100000 500000 1000000)
DOWNSAMPLE=(0 1)

OUTPUT_DIR="${1}"
[[ -d $OUTPUT_DIR ]] || (echo 'Invalid output dir' && exit 1)
HIC_SAMPLES=${@:2}
# All combos of hicrep params
for max_window_size in ${MAX_WINDOW_SIZE[@]}; do
for resolution in ${RESOLUTIONS[@]}; do
for h in ${HS[@]}; do
for downsample in ${DOWNSAMPLE[@]}; do
    # For each possible valid combination of parameters
    if [[ $downsample -eq 0 ]]; then
        downsample="" 
        downsample_str="original" 
    else
        downsample="--bDownSample"
        downsample_str="downsampled"
    fi
    [[ $max_window_size -lt $resolution ]] && continue
    # For all combinations of samples
    for sample_file_i in ${HIC_SAMPLES[@]}; do
    for sample_file_j in ${HIC_SAMPLES[@]}; do
        # sample 1 info
        sample_file_i="$(readlink -e ${sample_file_i})"
        sample_ID_i="$(basename $sample_file_i)"
        sample_ID_i="${sample_ID_i%%.mcool}"
        sample_name_i="$(echo $sample_ID_i | cut -d'.' -f1)"
        sample_params_i="$(echo $sample_ID_i | cut -d'.' -f2-)"
        # sample 2 info
        sample_file_j="$(readlink -e ${sample_file_j})"
        sample_ID_j="$(basename $sample_file_j)"
        sample_ID_j="${sample_ID_j%%.mcool}"
        sample_name_j="$(echo $sample_ID_j | cut -d'.' -f1)"
        sample_params_j="$(echo $sample_ID_j | cut -d'.' -f2-)"
        # Check if results are redundant or already exist
        comparison_dir="${downsample_str}/${max_window_size}/${h}/${resolution}"
        output_file="${OUTPUT_DIR}/${comparison_dir}/${sample_ID_i}.vs.${sample_ID_j}-hicrep.txt"
        redundant_file="${OUTPUT_DIR}/${comparison_dir}/${sample_ID_j}.vs.${sample_ID_i}-hicrep.txt"
        [[ -e ${output_file} || -e ${redundant_file} ]] && continue
        # Run hicrep
        mkdir -p "${OUTPUT_DIR}/${comparison_dir}"
        echo $output_file
        hicrep $downsample \
            --h $h \
            --dBPMax $max_window_size \
            --binSize $resolution \
            $sample_file_i \
            $sample_file_j \
            $output_file
    done
    done
done
done
done
done
