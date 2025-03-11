#!/bin/bash
# Mitochondrial Genome Assembly Pipeline v2024.12.1
# Author: WJJ

# ============================
# Usage Guide
# ============================
# Usage:
# Mitochondria_Assemble <fastq_1> <fastq_2> <species_name> <steps>
#
# Parameters:
#   <fastq_1>       : First paired-end sequencing read file (e.g., R1.fastq)
#   <fastq_2>       : Second paired-end sequencing read file (e.g., R2.fastq)
#   <species_name>  : Species name for directory and logging
#   <steps>         : Steps to execute (1=Mitoz, 2=NOVOPlasty, e.g., "12")
#
# Features:
# - Detailed time tracking for each step
# - Structured logging with timestamps
# - Error checking and validation
# ============================

# Enable strict error handling
set -euo pipefail
IFS=$'\n\t'

# Validate arguments
if [ "$#" -ne 4 ]; then
    echo "Error: Required 4 arguments. Usage:"
    echo "  $0 <fastq_1> <fastq_2> <species> <steps>"
    exit 1
fi

# Initialize variables
FASTQ1="$1"
FASTQ2="$2"
SPECIES_NAME="$3"
STEPS="$4"
BASE_DIR=$(pwd)
START_TIME=$(date +%s)

# Create directory structure
LOG_DIR="${BASE_DIR}/0-log"
SPECIES_DIR="${BASE_DIR}/${SPECIES_NAME}"
LOG_FILE="${LOG_DIR}/${SPECIES_NAME}_${STEPS}_assemble.log"

mkdir -p "$LOG_DIR" "$SPECIES_DIR"
echo "species_name=${SPECIES_NAME}" > "$LOG_FILE"

# Timing functions
get_timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}

format_duration() {
    local seconds=$1
    printf "%dh %dm %ds" $((seconds/3600)) $(( (seconds%3600)/60 )) $((seconds%60))
}

# Step execution tracking
declare -A STEP_TIMING=()

# MitoZ Assembly
run_mitoz() {
    local step_start=$(date +%s)
    echo -e "\n# MitoZ [$(get_timestamp)]" >> "$LOG_FILE"
    
    local step_dir="${SPECIES_DIR}/1-${SPECIES_NAME}_Mitoz"
    mkdir -p "${step_dir}"
    cd "${step_dir}"

    echo "Starting MitoZ assembly..." | tee -a "$LOG_FILE"
    
    mitoz all \
        --thread_number 20 \
        --clade "Arthropoda" \
        --genetic_code 5 \
        --species_name "${SPECIES_NAME}" \
        --fq1 "${FASTQ1}" \
        --fq2 "${FASTQ2}" \
        --insert_size 350 \
        --fastq_read_length 150 \
        --skip_filter \
        --data_size_for_mt_assembly 0,8 \
        --assembler megahit \
        --kmers 33 45 55 77 99 119 \
        --memory 40 \
        --slow_search \
        --requiring_taxa "Arthropoda" >> "$LOG_FILE"

    local step_end=$(date +%s)
    STEP_TIMING[1]=$((step_end - step_start))
    echo "MitoZ completed in $(format_duration ${STEP_TIMING[1]})" | tee -a "$LOG_FILE"
    cd "$BASE_DIR"
}

# NOVOPlasty Assembly
run_novoplasty() {
    local step_start=$(date +%s)
    echo -e "\n# NOVOPlasty [$(get_timestamp)]" >> "$LOG_FILE"
    
    local step_dir="${SPECIES_DIR}/2-${SPECIES_NAME}_NOVOPlasty"
    local config_dir="${BASE_DIR}/2-configs"
    local config_file="${config_dir}/NOVOPlasty_config_${SPECIES_NAME}.txt"
    
    mkdir -p "${step_dir}" "${config_dir}"
    cd "${step_dir}"

    # Generate config file
    cat > "$config_file" << EOF
Project:
-----------------------
Project name          = ${SPECIES_NAME}
Type                  = mito
Genome Range          = 12000-20000
K-mer                 = 33
Max memory            = 40
Extended log          = 0
Save assembled reads  = no
Seed Input            = ref.fasta
Extend seed directly  = no
Reference sequence    = ref.fasta
Variance detection    = no

Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 350
Platform              = illumina
Single/Paired         = PE
Forward reads         = ${FASTQ1}
Reverse reads         = ${FASTQ2}
EOF

    echo "Starting NOVOPlasty assembly..." | tee -a "$LOG_FILE"
    NOVOPlasty.pl -c "$config_file" >> "$LOG_FILE"
    rm "$config_file"

    local step_end=$(date +%s)
    STEP_TIMING[2]=$((step_end - step_start))
    echo "NOVOPlasty completed in $(format_duration ${STEP_TIMING[2]})" | tee -a "$LOG_FILE"
    cd "$BASE_DIR"
}

# Main execution
IFS='' read -r -a steps_array <<< $(echo "$STEPS" | grep -o .)
for step in "${steps_array[@]}"; do
    case "$step" in
        1) run_mitoz ;;
        2) run_novoplasty ;;
        *) echo "Warning: Skipping invalid step $step" | tee -a "$LOG_FILE" ;;
    esac
done

# Final report
TOTAL_TIME=$(( $(date +%s) - START_TIME ))
echo -e "\n# Execution Summary [$(get_timestamp)]" >> "$LOG_FILE"
for step in "${!STEP_TIMING[@]}"; do
    echo "Step $step: $(format_duration ${STEP_TIMING[$step]})" >> "$LOG_FILE"
done
echo "Total time: $(format_duration $TOTAL_TIME)" | tee -a "$LOG_FILE"

echo "Assembly completed. Log: ${LOG_FILE}"
