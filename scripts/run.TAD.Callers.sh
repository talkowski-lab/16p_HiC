# #!/bin/bash
# set -euo pipefail

###################################################
# TAD Calling params
###################################################
# hiTAD params
declare -rA HITAD_WEIGHTS=(["ICE"]="weight" ["Raw"]="RAW")
# cooltools insulation params
COOLTOOLS_WINDOW_SIZES="20 60 100" # numer of bins, not bp 
COOLTOOLS_MFVP=(0.33 0.66 0.9)
# COOLTOOLS_MFVP=( 0.66 0.9)
COOLTOOLS_THRESHOLD=(Li 0)
# COOLTOOLS_THRESHOLD=(Li)
declare -rA COOLTOOLS_WEIGHTS=(["ICE"]="weight" ["Raw"]="")

###################################################
# Functions
###################################################
help() {
    echo "USAGE: $(basename ${0}) [OPTIONS] {METHOD} sample1.mcool sample{2..N}.mcool
        -a | --anaconda-dir        # where conda is installed
        -r | --resolution          # resolutions at which to annotate TADs
        -o | --output-dir          # root results dir
        -h | --help                # print this message
            -l | --log-dir
            -p | --partition
            -m | --mem
            -t | --ntasks-per-node
            -c | --cpus
"
    exit 0
}

activate_conda() {
    # activate conda env with specific tools for each task
    case "$1" in
        hiTAD)     env_name="TADLib" ;;
        cooltools) env_name="cooltools" ;;
        *)         echo "Invalid conda env: $1" && exit 1 ;;
    esac
    echo "source ${CONDA_DIR}/etc/profile.d/conda.sh; conda activate ${env_name}"
}

get_sample_ID() {
    sample_ID="$(basename "${sample_file}")"
    sample_ID="${sample_ID%%.mcool}"
    sample_ID="${sample_ID%%.hic}"
    echo "${sample_ID}"
}

run_cmd() {
    cmd="${1}"
    job_name="${2}"
    case ${EVAL_METHOD} in 
        slurm)   
            sbatch                                      \
                --partition "${QUEUE}"                  \
                --ntasks="${NTASKS_PER_NODE}"           \
                --cpus-per-task="${CPUS}"               \
                --mem="${MEM_GB}G"                      \
                --job-name "${job_name}"                \
                --output   "${LOG_DIR}/${job_name}.out" \
                --error    "${LOG_DIR}/${job_name}.err" \
                --wrap="${cmd}"
            ;;
        txt)     
            echo "${CONDA_ENV_CMD}; ${cmd}" 
            ;;
        inplace) 
            time eval "${cmd}" 
            ;;
        *)       
            echo "Invalid launch method ${EVAL_METHOD}" && help && exit 1 
            ;;
    esac
}

run_hitad() {
    # input args
    method='hiTAD'
    mcool_file="$1"
    resolution="$2"
    # get sample ID
    mcool_file="$(readlink -e "${mcool_file}")"
    sample_ID="$(get_sample_ID "${mcool_file}")"
    # echo "${sample_ID}"
    # For all tool-specific param combos call TADs
    for weight_name in "${!HITAD_WEIGHTS[@]}"; do
        # log info
        param_dir="method_${method}/resolution_${resolution}/weight_${weight_name}"
        output_dir="${OUTPUT_DIR}/${param_dir}"
        output_TAD_file="${output_dir}/${sample_ID}-TAD.tsv"
        # Dont submit job if results already exist
        if [[ -a "${output_TAD_file}" ]]; then
            continue 
        fi
        # compose command + args to annotate TADs
        mkdir -p "${output_dir}"
        cmd="domaincaller --cpu-core ${THREADS} --weight-col ${HITAD_WEIGHTS[${weight_name}]} --logFile ${output_dir}/${sample_ID}-hitad.log --DI-output ${output_dir}/${sample_ID}-DI.tsv --output ${output_TAD_file} --uri ${mcool_file}::resolutions/${resolution}"
        # domaincaller
        #     --cpu-core   ${THREADS}
        #     --weight-col ${HITAD_WEIGHTS[${weight_name}]}
        #     --logFile    ${output_dir}/${sample_ID}-hitad.log
        #     --DI-output  ${output_dir}/${sample_ID}-DI.tsv
        #     --output     ${output_TAD_file}
        #     --uri ${mcool_file}::resolutions/${resolution}
        # run command
        job_name="$(echo "${param_dir}" | sed -e 's/\//./g').${sample_ID}"
        run_cmd "${cmd}" "${job_name}"
    done
}

run_cooltools_insulation() {
    # input args
    method='cooltools'
    mcool_file="$1"
    resolution="$2"
    # get sample ID
    mcool_file="$(readlink -e "${mcool_file}")"
    sample_ID="$(get_sample_ID "${mcool_file}")"
    # For all tool-specific param combos call TADs on this sample at this resol
    for weight_name in ${!COOLTOOLS_WEIGHTS[@]}; do
    for threshold in ${COOLTOOLS_THRESHOLD[@]}; do
    for mfvp in ${COOLTOOLS_MFVP[@]}; do
        # log info
        param_dir="method_${method}/resolution_${resolution}/weight_${weight_name}/threshold_${threshold}/mfvp_${mfvp}"
        # output files
        output_dir="${OUTPUT_DIR}/${param_dir}"
        output_TAD_file="${output_dir}/${sample_ID}-TAD.tsv"
        # Dont submit job if results already exist
        if [[ -a ${output_TAD_file} ]]; then
            continue 
        fi
        # Specify TAD calling cmd + options
        if [ -n "${COOLTOOLS_WEIGHTS[${weight_name}]}" ]; then
            weight_flag="--clr-weight-name ${COOLTOOLS_WEIGHTS[${weight_name}]} "
        else
            weight_flag=""
        fi
        mkdir -p "${output_dir}"
        cmd="cooltools insulation ${weight_flag}--verbose --nproc ${THREADS} --append-raw-scores --window-pixels --min-frac-valid-pixels ${mfvp} --ignore-diags 2 --threshold ${threshold} --output ${output_TAD_file} ${mcool_file}::resolutions/${resolution} ${COOLTOOLS_WINDOW_SIZES[*]}"
        # cooltools insulation 
        #     ${weight_flag}
        #     --verbose
        #     --nproc ${THREADS}
        #     --append-raw-scores                                   
        #     --window-pixels                                       
        #     --min-frac-valid-pixels ${mfvp}                       
        #     --ignore-diags 2 
        #     --threshold ${threshold}                              
        #     --output ${output_TAD_file} 
        #     ${mcool_file}::resolutions/${resolution} 
        #     ${COOLTOOLS_WINDOW_SIZES[*]}"
        # Launch as slurm job if specified
        job_name="$(echo "${param_dir}" | sed -e 's/\//./g').${sample_ID}"
        run_cmd "${cmd}" "${job_name}"
    done
    done
    done
}

main() {
    # Activate conda env now if running locally (not submitting jobs)
    [[ ${EVAL_METHOD} == 'inplace' ]] && eval "${CONDA_ENV_CMD}"
    # Loop over all samples at all resolutions and call tads
    for resolution in ${RESOLUTIONS[@]}; do
        for sample_file in ${HIC_SAMPLES[@]}; do
            sample_ID="$(get_sample_ID "${sample_file}")"
            echo "${sample_ID}"
            case ${METHOD} in
                hiTAD)     run_hitad "${sample_file}" "${resolution}" ;;
                cooltools) run_cooltools_insulation "${sample_file}" "${resolution}" ;;
                *)         echo "Invalid method: ${METHOD}" && exit 1 ;;
            esac
        done
    done
}

###################################################
# Handle CLI args
###################################################
# Default script arguments
# BASE_DIR="/data/talkowski/Samples/16p_HiC"
BASE_DIR="./"
OUTPUT_DIR="${BASE_DIR}/results/TADs"
RESOLUTIONS=(100000 50000 25000 10000)
# Default SLURM params
USE_SLURM=0
PARTITION="normal"
MEM_GB=30
NTASKS_PER_NODE=1
CPUS=2
THREADS=2
CONDA_DIR="$(conda info --base)"
# Handle CLI args
[[ $? -ne 0 ]] && echo "No Args" && exit 1
VALID_ARGS=$(getopt -o ho:l:r:n:c:p:m:a:t: --long help,output-dir,log-dir,resolution,nstasks-per-node,cpus,partition,mem,anaconda-dir,threads -- "$@")
eval set -- "${VALID_ARGS}"
while [ : ]; do
    echo "${1} ||| ${2}"
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
        -t|--threads)
            THREADS="${2}" 
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
        -n|--ntasks-per-node)
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

###################################################
# Main 
###################################################
METHOD="${1}"
HIC_SAMPLES=${@:2}
echo "
Using TAD caller:       ${METHOD}
Using resolution(s):    ${RESOLUTIONS[*]}
Using output directory: ${OUTPUT_DIR}
Threads:                ${THREADS}"
CONDA_ENV_CMD="$(activate_conda "${METHOD}")"
OUTPUT_DIR="$(readlink -e "${OUTPUT_DIR}")"
mkdir -p "${OUTPUT_DIR}"
main 
