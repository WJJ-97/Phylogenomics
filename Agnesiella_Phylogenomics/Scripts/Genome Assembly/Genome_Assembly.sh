#!/bin/bash
# Genome Assembly Pipeline v2024.12.1
# Author: WJJ

# ============================
# Global Configuration
# ============================
set -euo pipefail
IFS=$'\n\t'

# ============================
# Function Definitions
# ============================

# Unified logging with timestamp and duration tracking
log() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] $1" | tee -a "$LOG_FILE"
}

# Track execution time for pipeline steps
track_step() {
    local step_name=$1
    local start_time=$2
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log "STEP COMPLETED: $step_name | Duration: ${duration}s"
}

# Safe directory change with error handling
safe_cd() {
    cd "$1" || {
        log "FATAL: Failed to enter directory: $1"
        exit 1
    }
}

# Dependency verification
check_dependency() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "FATAL: Missing required dependency: $1"
        exit 1
    }
}

# Cleanup temporary files
clean_temp_files() {
    local step=$1
    log "Cleaning temporary files for step $step"
    find "$BASE_PATH" -name '*tmp*' -delete 2>/dev/null || true
}

# Validate input parameters
validate_parameters() {
    [[ "$#" -ge 3 ]] || {
        echo "Usage: $0 <forward.fq> <reverse.fq> <species> [steps] [threads] [memory]"
        echo "Example: $0 read1.fq read2.fq Drosophila 12345678 16 64G"
        exit 1
    }

    [[ -f "$1" && -f "$2" ]] || {
        log "FATAL: Input files not found: $1 and $2"
        exit 1
    }

    STEPS=${4:-12345678}
    [[ "$STEPS" =~ ^[1-8]+$ ]] || {
        log "FATAL: Invalid step format: $STEPS (Must be 1-8 combinations)"
        exit 1
    }

    # Deduplicate and sort steps
    STEPS=$(echo "$STEPS" | fold -w1 | sort -u | tr -d '\n')
}

# ============================
# Pipeline Modules
# ============================

run_quality_trim() {
    local start=$(date +%s)
    log "STEP 1: QUALITY TRIMMING - STARTED"
    
    mkdir -p "${BASE_PATH}/1-${SPECIES_NAME}_fastp"
    safe_cd "${BASE_PATH}/1-${SPECIES_NAME}_fastp"

    fastp --in1 "$FORWARD_READS" \
          --in2 "$REVERSE_READS" \
          --out1 "${PREFIX1}.fastp.fq.gz" \
          --out2 "${PREFIX2}.fastp.fq.gz" \
          --html "${SPECIES_NAME}_report.html" \
          --thread "$THREADS" \
          --detect_adapter_for_pe \
          --trim_front1 1 --trim_front2 1 \
          --trim_tail1 1 --trim_tail2 1 \
          --length_required 20 \
          2>> "$LOG_FILE"

    clean_temp_files 1
    track_step "Quality Trimming" $start
}

run_normalization() {
    local start=$(date +%s)
    log "STEP 2: READ NORMALIZATION - STARTED"
    
    mkdir -p "${BASE_PATH}/2-${SPECIES_NAME}_normalize"
    safe_cd "${BASE_PATH}/2-${SPECIES_NAME}_normalize"

    bbnorm.sh in1="../1-${SPECIES_NAME}_fastp/${PREFIX1}.fastp.fq.gz" \
              in2="../1-${SPECIES_NAME}_fastp/${PREFIX2}.fastp.fq.gz" \
              out1="${PREFIX1}.nor.fq.gz" \
              out2="${PREFIX2}.nor.fq.gz" \
              target=10 min=2 \
              threads="$THREADS" -Xmx"${MEMORY}" \
              2>> "$LOG_FILE"

    clean_temp_files 2
    track_step "Read Normalization" $start
}

run_spades_assembly() {
    local start=$(date +%s)
    log "STEP 3: SPADES ASSEMBLY - STARTED"
    
    mkdir -p "${BASE_PATH}/3-${SPECIES_NAME}_spades"
    safe_cd "${BASE_PATH}/3-${SPECIES_NAME}_spades"

    spades.py -1 "../2-${SPECIES_NAME}_normalize/${PREFIX1}.nor.fq.gz" \
             -2 "../2-${SPECIES_NAME}_normalize/${PREFIX2}.nor.fq.gz" \
             -o "assembly" \
             -k 21,33,55,77 \
             -t "$THREADS" \
             2>> "$LOG_FILE"

    clean_temp_files 3
    track_step "SPAdes Assembly" $start
}

# ============================
# Main Execution Flow
# ============================

main() {
    # Initialize execution tracking
    local start_time=$(date +%s)
    trap 'log "WARNING: User interrupted execution"; exit 1' SIGINT

    # Parameter handling
    validate_parameters "$@"
    FORWARD_READS=$(realpath "$1")
    REVERSE_READS=$(realpath "$2")
    SPECIES_NAME="$3"
    THREADS=${5:-16}
    MEMORY=${6:-64G}
    PREFIX1=$(basename "${FORWARD_READS%%.*}")
    PREFIX2=$(basename "${REVERSE_READS%%.*}")

    # System configuration
    BASE_PATH=$(realpath "./${SPECIES_NAME}")
    LOG_DIR="${BASE_PATH}/0-log"
    mkdir -p "$LOG_DIR"
    LOG_FILE="${LOG_DIR}/${SPECIES_NAME}_assembly.log"

    # Dependency verification
    check_dependency fastp
    check_dependency bbnorm.sh
    check_dependency spades.py

    # Execute pipeline steps
    for step in $(echo "$STEPS" | fold -w1); do
        case $step in
            1) run_quality_trim ;;
            2) run_normalization ;;
            3) run_spades_assembly ;;
            4) run_heterozygosity_reduction ;;
            5) run_mapping ;;
            6) run_scaffolding ;;
            7) run_gap_closing ;;
            8) run_busco ;;
            *) log "WARNING: Ignored invalid step $step" ;;
        esac
    done

    # Finalization
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log "PIPELINE COMPLETED: Total runtime $((duration / 3600))h $(( (duration % 3600) / 60 ))m $((duration % 60))s"
}

# ============================
# Entry Point
# ============================
main "$@"
