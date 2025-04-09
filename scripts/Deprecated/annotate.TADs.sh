#!/bin/bash
set -euo pipefail
# Locations
BASE_DIR="$(pwd)"
DISTILLER_FILE="${BASE_DIR}/distiller-nf/distiller.nf"
JUICER_JAR="/data/talkowski/Samples/WAPL_NIPBL/HiC/distiller-nf/juicer_tools_1.22.01.jar"
JAVA_CMD="java -jar -Xmx48000m -Djava.awt.headless=true -jar ${JUICER_JAR}"
# arrowhead params
METHODS=(arrowhead)
RESOLUTIONS=(100000 50000 10000 5000)
WINDOW_SIZES=(4000 2000)
BALANCINGS=(KR VC NONE)
WEIGHTS=(weight)
WEIGHT_NAMES=(ICE)
# SLURM params
QUEUE="short"
THREADS=2
MEM_GB=30
# Initialize loop
jobs_done=0
jobs_skipped=0
OUTPUT_DIR="${1}"
[[ -d $OUTPUT_DIR ]] || (echo 'Invalid output dir' && exit 1)
LOG_DIR="${2}"
[[ -d $LOG_DIR ]] || (echo 'Invalid log dir' && exit 1)
HIC_SAMPLES=${@:3}
# Run a slurm job calling TADs for each samples
# example ./run.arrowhead.sh ./samples/HiC.sample{1..6}.hic 
for sample_file in ${HIC_SAMPLES[@]}; do
    sample_file="$(readlink -e ${sample_file})"
    sample_ID="$(basename $sample_file)"
    sample_ID="${sample_ID%%.hic}"
    sample_ID="${sample_ID%%.mcool}"
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
echo "Ran ${job_num} jobs, skipped ${jobs_skipped}"
echo "${job_num} jobs submitted on queue ${QUEUE} with ${MEM_GB}Gb, ${THREADS} threads per job"

run_hitad() {
    domaincaller                      \
        --logFile ${log_file}         \
        --cpu-core $THREADS           \
        --weight-col ${weight}        \
        --output  ${output_TAD_file}  \
        --DI-output ${output_DI_file} \
        --uri "${sample_file}::resolutions/${resolution}"
}
run_arrowhead() {
    echo None
}
run_cooltools_insulation() {
    cooltools insulation \
        --verbose
}
run_slurm_cmd() {
    sbatch                                      \
        --partition ${QUEUE}                    \
        --nodes=${NODES}                        \
        --ntasks=${THREADS}                     \
        --mem="${MEM_GB}G"                      \
        --job-name "${job_name}"                \
        --output   "${LOG_DIR}/${job_name}.out" \
        --error    "${LOG_DIR}/${job_name}.err" \
        --wrap="${cmd}"
}
#!/bin/bash
OUTPUT_DIR="${1}"
LOG_DIR="${2}"
HIC_SAMPLES=${@:3}
# loop over all param combos
for resolution in ${RESOLUTIONS[@]}; do
for sample_file in ${HIC_SAMPLES[@]}; do
    case $METHOD in
        hiTAD)
            ;;
        arrowhead)
            ;;
        cooltools)
            ;;
        *)
            ;;
    esac
for i in ${!WEIGHTS[@]}; do
    weight=${WEIGHTS[i]}
    weight_name=${WEIGHT_NAMES[i]}
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
done
done
done
done
