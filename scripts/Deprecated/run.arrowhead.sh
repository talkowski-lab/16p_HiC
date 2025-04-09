#!/bin/bash
set -euo pipefail
# Locations
# DISTILLER_DIR="/data/talkowski/Samples/WAPL_NIPBL/HiC/distiller-nf"
# LOG_DIR="${DISTILLER_DIR}/logs"
JUICER_JAR="/data/talkowski/Samples/WAPL_NIPBL/HiC/distiller-nf/juicer_tools_1.22.01.jar"
JAVA_CMD="java -jar -Xmx48000m -Djava.awt.headless=true -jar ${JUICER_JAR}"
# arrowhead params
RESOLUTIONS=(100000 50000 10000 5000)
WINDOW_SIZES=(4000 2000)
BALANCINGS=(KR VC NONE)
# RESOLUTIONS=(100000)
# WINDOW_SIZES=(2000)
# BALANCINGS=(KR)
METHODS=(arrowhead)
# SLURM params
QUEUE="short"
THREADS=2
MEM_GB=30
# Input files
i=0
OUTPUT_DIR="${1}"
LOG_DIR="${2}"
HIC_SAMPLES=${@:3}
# example ./run.arrowhead.sh ./samples/HiC.sample{1..6}.hic 
for sample_file in ${HIC_SAMPLES[@]}; do
    sample_file="$(readlink -e ${sample_file})"
    sample_ID="$(basename $sample_file)"
    sample_ID="${sample_ID%%.hic}"
    # echo $sample_ID
    for resolution in ${RESOLUTIONS[@]}; do
    for window_size in ${WINDOW_SIZES[@]}; do
    for balancing in ${BALANCINGS[@]}; do
    for method in ${METHODS[@]}; do
        job_name="${method}.${balancing}.${window_size}.${sample_ID}.${resolution}"
        log_file="${LOG_DIR}/${job_name}.log"
        # Dont submit job if results already exist
        output_dir="${OUTPUT_DIR}/${method}/${balancing}/${window_size}/${sample_ID}"
        output_file="${output_dir}/${resolution}_blocks.bedpe"
        [[ -s $output_file ]] && continue || i=$(( i+1 ))
        # continue
        mkdir -p "${output_dir}"
        basename "${log_file}"
        echo $sample_file
        # continue
        # allow multiple methods
        case "${method}" in
            arrowhead) 
                cmd="source /PHShome/sr1068/miniforge3/etc/profile.d/conda.sh
conda activate dist2
${JAVA_CMD} arrowhead    \
    --threads ${THREADS} \
    -k ${balancing}      \
    -m ${window_size}    \
    -r ${resolution}     \
       ${sample_file}    \
       ${output_file}"
                ;;
            *)
                echo "Error: ${method} not implemented"
                exit 1
                ;;
        esac
        # Run job on slurm with specified params
        sbatch                         \
            --job-name  "${job_name}" \
            --output    "${log_file}" \
            --partition "${QUEUE}"    \
            --ntasks    "${THREADS}"  \
            --mem       "${MEM_GB}G"  \
            --wrap="${cmd}"
        # $cmd
    done
    done
    done
    done
done
echo "$i jobs submitted on ${QUEUE} queue, ${MEM_GB}Gb. ${THREADS} threads per job"
