#!/bin/bash
# Args are hic/mcool files to convert
# SLURM params
QUEUE="short"
THREADS=12
MEM_GB=40
# Inputs 
OUTPUT_DIR="${1}"
[[ -d $OUTPUT_DIR ]] || (echo 'Invalid output dir' && exit 1)
LOG_DIR="${2}"
[[ -d $LOG_DIR ]] || (echo 'Invalid log dir' && exit 1)
HIC_SAMPLES=${@:3}
i=0
for sample_file in ${HIC_SAMPLES[@]}; do
    filename="$(basename ${sample_file})"
    sample_ID="${filename%.*}"
    extension="${filename##*.}"
    case $extension in 
        mcool)
            output_file="${OUTPUT_DIR}/${sample_ID}.hic"
            ;;
        hic)
            output_file="${OUTPUT_DIR}/${sample_ID}.mcool"
            ;;
        *)
            echo "Invalid file extension: $extension"
            continue
            ;;
    esac
    job_name="${sample_ID}-convert_format"
    log_file="${LOG_DIR}/${job_name}.log"
    [[ -s $output_file ]] && continue || i=$(( i+1 ))
    echo $sample_ID
    # Run job on slurm with specified params
    sbatch                        \
        --job-name  "${job_name}" \
        --output    "${log_file}" \
        --partition "${QUEUE}"    \
        --ntasks    "${THREADS}"  \
        --mem       "${MEM_GB}G"  \
        --wrap="hictk convert --threads $THREADS $(readlink -e ${sample_file}) ${output_file}"
done
echo "$i jobs submitted on ${QUEUE} queue, ${MEM_GB}Gb. ${THREADS} threads per job"
