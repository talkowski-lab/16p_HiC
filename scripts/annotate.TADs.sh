# #!/bin/bash
set -euo pipefail
# Arrowhead params
JUICER_JAR="/data/talkowski/Samples/WAPL_NIPBL/HiC/distiller-nf/juicer_tools_1.22.01.jar"
JAVA_CMD="java -jar -Xmx48000m -Djava.awt.headless=true -jar ${JUICER_JAR}"
ARROWHEAD_WINDOW_SIZES=(4000 2000)
ARROWHEAD_BALANCINGS=(KR VC NONE)
# hiTAD params
declare -rA HITAD_WEIGHTS=(["ICE"]="weight"  ["Raw"]="RAW")
# cooltools insulation params
COOLTOOLS_WINDOW_SIZES="20 60 100" # numer of bins, not bp 
COOLTOOLS_MFVP=(0.33 0.66 0.9)
COOLTOOLS_THRESHOLD=(Li 0)
declare -rA COOLTOOLS_WEIGHTS=(["ICE"]="weight"  ["Raw"]="")
# Functions
help() {
    echo "USAGE: $(basename $0) [OPTIONS] {METHOD} sample1.mcool sample{2..N}.mcool 
            -s | --use-slurm
            -r | --resolution
            -o | --output-dir
            -l | --log-dir
            -h | --help"
    exit 0
}
activate_conda() {
    # activate conda env with specific tools for each task
    case "$1" in
        arrowhead) env_name="dist2" ;;
        cooltools) env_name="cooltools" ;;
        hiTAD)     env_name="TADLib" ;;
        *)         echo "Invalid conda env: $1" && exit 1 ;;
    esac
    echo "source ${CONDA_DIR}/etc/profile.d/conda.sh; conda activate ${env_name}"
}
get_sample_ID() {
    sample_ID="$(basename ${sample_file})"
    sample_ID="${sample_ID%%.mcool}"
    sample_ID="${sample_ID%%.hic}"
    echo "${sample_ID}"
}
run_slurm_cmd() {
    cmd="$1"
    job_name="$2"
    sbatch                                      \
        --partition ${PARTITION}                \
        --ntasks=${NTASKS_PER_NODE}             \
        --cpus-per-task=${CPUS}                 \
        --mem=${MEM_GB}G                        \
        --job-name "${job_name}"                \
        --output   "${LOG_DIR}/${job_name}.out" \
        --error    "${LOG_DIR}/${job_name}.err" \
        --wrap="${cmd}"
}
run_arrowhead() {
    hic_file="$1"
    output_file="$2"
    balancing="$3"
    resolution="$4"
    window_size="$5"
    cmd="${CONDA_ENV_CMD}
    ${JAVA_CMD} arrowhead
        --threads ${NTASKS_PER_NODE}
        -k ${balancing}
        -r ${resolution}
        -m ${window_size}
           ${hic_file}
           ${output_file}"
}
run_hitad() {
    # input args
    method='hiTAD'
    mcool_file="$1"
    resolution="$2"
    # keep track of jobs submitted
    jobs_num=0
    jobs_done=0
    jobs_skipped=0
    # get sample ID
    mcool_file="$(readlink -e ${mcool_file})"
    sample_ID="$(get_sample_ID "${mcool_file}")"
    # For all tool-specific param combos call TADs
    for weight_name in ${!HITAD_WEIGHTS[@]}; do
        jobs_num=$(( jobs_num+1 ))
        # log info
        param_dir="${method}/${resolution}/${weight_name}"
        output_dir="${OUTPUT_DIR}/${param_dir}"
        output_TAD_file="${output_dir}/${sample_ID}-TAD.tsv"
        # Dont submit job if results already exist
        if [[ -a $output_TAD_file ]]; then
            jobs_skipped=$(( jobs_skipped+1 )) 
            continue 
        fi
        # Run command
        mkdir -p "${output_dir}"
        cmd="domaincaller                                     \
            --cpu-core   ${NTASKS_PER_NODE}                           \
            --weight-col ${HITAD_WEIGHTS[$weight_name]}       \
            --logFile    ${output_dir}/${sample_ID}-hitad.log \
            --DI-output  ${output_dir}/${sample_ID}-DI.tsv    \
            --output     ${output_TAD_file}                   \
            --uri ${mcool_file}::resolutions/${resolution}"
        # make sure conda env is activated in the job if using slurm
        # echo "${cmd}"
        if [[ $USE_SLURM -eq 1 ]]; then
            job_name="$(echo $param_dir | sed -e 's/\//./g').${sample_ID}"
            run_slurm_cmd "${CONDA_ENV_CMD}; ${cmd}" "${job_name}"
            jobs_done=$(( jobs_done+1 ))
        else 
            eval $cmd 
            [[ -a $output_TAD_file ]] && jobs_done=$(( jobs_done+1 ))
        fi
    done
    echo -e "${sample_ID}\t${jobs_num}\t${jobs_done}\t${jobs_skipped}"
}
run_cooltools_insulation() {
    # input args
    method='cooltools'
    mcool_file="$1"
    resolution="$2"
    # keep track of jobs submitted
    jobs_total=0
    jobs_num=0
    jobs_done=0
    jobs_skipped=0
    # get sample ID
    mcool_file="$(readlink -e ${mcool_file})"
    sample_ID="$(get_sample_ID "${mcool_file}")"
    # For all tool-specific param combos call TADs on this sample at this resol
    for weight_name in ${!COOLTOOLS_WEIGHTS[@]}; do
    for threshold in ${COOLTOOLS_THRESHOLD[@]}; do
    for mfvp in ${COOLTOOLS_MFVP[@]}; do
        jobs_num=$(( jobs_num+1 ))
        # log info
        param_dir="${method}/${resolution}/${weight_name}/${threshold}/${mfvp}"
        # output files
        output_dir="${OUTPUT_DIR}/${param_dir}"
        output_TAD_file="${output_dir}/${sample_ID}-TAD.tsv"
        # Dont submit job if results already exist
        if [[ -a $output_TAD_file ]]; then
            jobs_skipped=$(( jobs_skipped+1 )) 
            continue 
        fi
        # Specify TAD calling cmd + options
        if [ -n "${COOLTOOLS_WEIGHTS[${weight_name}]}" ]; then
            weight_flag="--clr-weight-name ${COOLTOOLS_WEIGHTS[${weight_name}]} "
        else
            weight_flag=""
        fi
        mkdir -p "${output_dir}"
        cmd="cooltools insulation ${weight_flag}--verbose
                --nproc ${NTASKS_PER_NODE}                                   
                --append-raw-scores                                   
                --window-pixels                                       
                --min-frac-valid-pixels ${mfvp}                       
                --ignore-diags 2 
                --threshold ${threshold}                              
                --output ${output_TAD_file} 
                ${mcool_file}::resolutions/${resolution} 
                ${COOLTOOLS_WINDOW_SIZES[@]}"
        # Launch as slurm job if specified
        if [[ $USE_SLURM -eq 1 ]]; then
            job_name="$(echo $param_dir | sed -e 's/\//./g').${sample_ID}"
            run_slurm_cmd "${CONDA_ENV_CMD}; ${cmd}" "${job_name}"
        else 
            eval $cmd 
            [[ -a $output_TAD_file ]] && jobs_done=$(( jobs_done+1 ))
        fi
    done
    done
    done
    echo -e "${sample_ID}\t${jobs_num}\t${jobs_done}\t${jobs_skipped}"
}
main() {
    # Activate conda env now if running locally (not submitting jobs)
    [[ $USE_SLURM -eq 0 ]] && eval "${CONDA_ENV_CMD}"
    # Loop over all samples at all resolutions and call tads
    echo -e "SampleID\tJobs Submitted\tJobs Skipped"
    for sample_file in ${HIC_SAMPLES[@]}; do
    for resolution in ${RESOLUTIONS[@]}; do
        sample_ID="$(get_sample_ID "${sample_file}")"
        case $METHOD in
            arrowhead)
                run_arrowhead ${sample_file} ${resolution}
                ;;
            hiTAD)
                run_hitad ${sample_file} ${resolution}
                ;;
            cooltools)
                run_cooltools_insulation ${sample_file} ${resolution}
                ;;
            *)
                echo "Invalid method: $METHOD" && exit 1
                ;;
        esac
    done
    done
    # echo "All jobs submitted on queue ${PARTITION} with ${MEM_GB}Gb, ${NTASKS_PER_NODE} threads per job"
}
# Default script arguments
BASE_DIR="/data/talkowski/Samples/16p_HiC"
OUTPUT_DIR="${BASE_DIR}/results/TADs"
LOG_DIR="${BASE_DIR}/slurm.logs"
RESOLUTIONS=(100000 50000 25000 10000)
# Default SLURM params
USE_SLURM=0
PARTITION="short"
MEM_GB=30
NTASKS_PER_NODE=2
CPUS=2
# Handle CLI args
[[ $? -ne 0 ]] && echo "No Args" && exit 1
VALID_ARGS=$(getopt -o ho:l:r:t:c:p:m:a: --long help,output-dir,log-dir,resolution,nstasks-per-node,cpus,partition,mem,anaconda-dir -- "$@")
eval set -- "$VALID_ARGS"
while [ : ]; do
    case "$1" in
        -a|--anaconda-dir)
            CONDA_DIR="${2}"
            shift 2
            ;;
        -r|--resolution)
            RESOLUTIONS=($(echo "${2}" | cut -d',' --output-delimiter=' ' -f1-))
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="${2}" 
            shift 2
            ;;
        -s|--use-slurm)
            USE_SLURM=1
            shift 
            ;;
        -m|--mem)
            MEM_GB="${2}" 
            shift 2
            ;;
        -t|--ntasks-per-node)
            NTASKS_PER_NODE="${2}" 
            shift 2
            ;;
        -c|--cpus)
            CPUS="${2}" 
            shift 2
            ;;
        -p|--partition)
            PARTITION="${2}" 
            shift 2
            ;;
        -l|--log-dir)
            LOG_DIR="${2}" 
            shift 2
            ;;
        -a|--anaconda-dir)
            CONDA_DIR="${2}"
            shift 2
            ;;
        -h|--help) 
            help 
            ;;
        --)
            shift 
            break
            ;;
    esac
done
# Print args
METHOD="${1}"
HIC_SAMPLES=${@:2}
echo "
Using TAD caller:       ${METHOD}
Using resolution(s):    ${RESOLUTIONS[@]}
Using output directory: ${OUTPUT_DIR}"
if [[ $USE_SLURM -eq 1 ]]; then
echo "Submitting each TAD annotation as a slurm job with the following params
Log directory:           ${LOG_DIR}
Submitting to partition: ${PARTITION}
Mem per job:             ${MEM_GB}Gb
NTasks per node per job: ${NTASKS_PER_NODE}
CPUs per task per job:   ${CPUS}"
else
    echo "Not using SLURM, running TAD annotations in current shell"
fi 
# Run 
CONDA_ENV_CMD="$(activate_conda ${METHOD})"
mkdir -p ${OUTPUT_DIR}
OUTPUT_DIR="$(readlink -e ${OUTPUT_DIR})"
main 
