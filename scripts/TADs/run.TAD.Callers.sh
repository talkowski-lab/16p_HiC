#!/bin/bash
# set -eo pipefail

###################################################
# TAD Calling params
###################################################
# hiTAD params
declare -rA HITAD_WEIGHTS=(["ICE"]="weight" ["Raw"]="RAW")
# cooltools insulation params
# COOLTOOLS_WINDOW_SIZES="20 60 100" # numer of bins, not bp 
# COOLTOOLS_MFVP=(0.33 0.66 0.9)
# COOLTOOLS_THRESHOLD=(Li 0)
COOLTOOLS_WINDOW_SIZES="60 100" # numer of bins, not bp 
COOLTOOLS_MFVP=(0.66 0.9)
COOLTOOLS_THRESHOLD=(Li)
declare -rA COOLTOOLS_WEIGHTS=(["ICE"]="weight" ["Raw"]="")

###################################################
# Functions
###################################################
help() {
    echo "USAGE: $(basename ${0}) [OPTIONS] {METHOD} sample{1..N}.mcool
    MEHTODS: cooltools, hiTAD
        [OPTIONS]
        -a  # where conda is installed
        -r  # resolutions at which to annotate TADs
        -o  # output dir
        -h  # print this message
        [OPTIONS] (Slurm)
        -l  # log-dir
        -p  # partition
        -m  # mem
        -t  # ntasks-per-node
        -c  # cpus
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
            echo "${cmd}" 
            # echo "${cmd}" >> ${CMD_TXT_FILE}
            ;;
        inplace) 
            time eval "${cmd}" 
            ;;
        *)       
            echo "Invalid launch method ${EVAL_METHOD}" && help
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
    mkdir -p "${OUTPUT_DIR}"
    # Activate conda env now if running locally (not submitting jobs)
    case ${EVAL_METHOD} in 
        slurm)   continue ;;
        txt) 
            CMD_TXT_FILE="${OUTPUT_DIR}/${METHOD}-TAD.calling.cmds.txt"
            rm "${CMD_TXT_FILE}"
            ;;
        inplace) 
            CONDA_ENV_CMD="$(activate_conda "${METHOD}")"
            eval "${CONDA_ENV_CMD}" 
            ;;
        *) echo "Invalid launch method ${EVAL_METHOD}" && help && exit 1 ;;
    esac
    # Loop over all samples at all resolutions and call tads
    for resolution in ${RESOLUTIONS[@]}; do
        for sample_file in ${HIC_SAMPLES[@]}; do
            sample_ID="$(get_sample_ID "${sample_file}")"
            echo "${sample_ID} @ ${resolution}"
            echo "================================"
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
EVAL_METHOD='txt'
RESOLUTIONS=(100000 50000 25000 10000)
# Default SLURM params
QUEUE="normal"; MEM_GB=30; NTASKS_PER_NODE=1; CPUS=2; THREADS=2; LOG_DIR="${BASE_DIR}/slurm.logs"
CONDA_DIR="$(conda info --base)"
[[ $# -eq 0 ]] && echo "No Args" && help
while getopts "ho:l:r:n:c:q:m:a:t:e:" flag; do
    case ${flag} in 
        h) help && exit 0 ;;
        e) EVAL_METHOD="${OPTARG}" ;;
        a) CONDA_DIR="${OPTARG}" ;;
        r) RESOLUTIONS=($(echo "${OPTARG}" | cut -d',' --output-delimiter=' ' -f1-)) ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        # SLURM params 
        t) THREADS="${OPTARG}" ;;
        m) MEM_GB="${OPTARG}" ;;
        n) NTASKS_PER_NODE="${OPTARG}" ;;
        c) CPUS="${OPTARG}" ;;
        q) QUEUE="${OPTARG}" ;;
        l) LOG_DIR="${OPTARG}" ;;
        *) echo "Invalid flag ${flag}" && help && exit 1 ;;
    esac
done
shift $(( OPTIND-1 ))

###################################################
# Main 
###################################################
METHOD="${1}"
HIC_SAMPLES=${@:2}
OUTPUT_DIR="$(readlink -e "${OUTPUT_DIR}")"
echo "
Using TAD caller:       ${METHOD}
Using resolution(s):    ${RESOLUTIONS[*]}
Using output directory: ${OUTPUT_DIR}
Threads:                ${THREADS}"
main 
