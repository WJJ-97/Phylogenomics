#!/bin/bash
# Genome Annotation Pipeline v2024.12.1
# Author: WJJ

# ============================
# Global Configuration
# ============================
set -euo pipefail
IFS=$'\n\t'

# ============================
# Function Definitions
# ============================

# Enhanced logging with execution time tracking
log() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] $1" | tee -a "$LOG_FILE"
}

# Time tracking wrapper for pipeline steps
track_execution() {
    local step_name=$1
    local start_time=$(date +%s)
    "$@"
    local exit_code=$?
    local duration=$(( $(date +%s) - start_time ))
    log "STEP COMPLETED: ${step_name} - Duration: $((duration / 3600))h $(( (duration % 3600) / 60 ))m $((duration % 60))s"
    return $exit_code
}

# Error handling for directory changes
safe_cd() {
    cd "$1" || {
        log "FATAL: Failed to enter directory: $1"
        exit 1
    }
}

# Dependency checker with version validation
check_dependency() {
    local cmd=$1
    local min_version=${2:-}
    command -v "$cmd" >/dev/null 2>&1 || {
        log "FATAL: Missing required dependency: $cmd"
        exit 1
    }
    
    [[ -n "$min_version" ]] && {
        local version=$($cmd --version 2>&1 | head -n1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?')
        [[ "$(printf "%s\n" "$min_version" "$version" | sort -V | head -n1)" == "$min_version" ]] || {
            log "FATAL: $cmd version $min_version+ required. Found $version"
            exit 1
        }
    }
}

# ============================
# Pipeline Steps
# ============================

run_repeat_masking() {
    log "Initiate repeat masking procedure"
    mkdir -p "${BASE_PATH}/2-mask"
    safe_cd "${BASE_PATH}/2-mask"

    windowmasker -mk_counts -in "$INPUT_GENOME" \
        -out "${SPECIES_NAME}_masking_library.ustat" \
        -mem "${MEMORY}" -sformat obinary 2>> "$LOG_FILE"

    windowmasker -ustat "${SPECIES_NAME}_masking_library.ustat" \
        -dust T -in "$INPUT_GENOME" \
        -out "${SPECIES_NAME}_genome.masked.fa" \
        -outfmt fasta 2>> "$LOG_FILE"

    clean_intermediates 1
}

# ... Similar modifications for other pipeline functions ...

# ============================
# Main Execution
# ============================

main() {
    local start_time=$(date +%s)
    trap 'log "Pipeline interrupted by user"; exit 1' SIGINT SIGTERM

    validate_parameters "$@"
    initialize_variables "$@"
    setup_environment
    
    # Execute pipeline steps with time tracking
    for (( i=0; i<${#STEPS}; i++ )); do
        case "${STEPS:i:1}" in
            1) track_execution "Repeat Masking" run_repeat_masking ;;
            2) track_execution "BRAKER Annotation" run_braker ;;
            3) track_execution "GeMoMa Annotation" run_gemoma ;;
            4) track_execution "TE Annotation" run_te_annotation ;;
        esac
    done

    finalize_execution "$start_time"
}

# ============================
# Helper Functions
# ============================

initialize_variables() {
    INPUT_GENOME=$(realpath "$1")
    SPECIES_NAME="$2"
    THREADS=${4:-8}
    MEMORY=${5:-40G}
    EMAIL=${6:-}
    BASE_PATH=$(realpath ./)
    LOG_DIR="${BASE_PATH}/0-log"
    LOG_FILE="${LOG_DIR}/${SPECIES_NAME}_${STEPS}.log"
}

setup_environment() {
    mkdir -p "$LOG_DIR"
    check_dependency windowmasker 2.2.28
    check_dependency braker.pl 2.1.6
    check_dependency java 11.0
    check_dependency earlGrey 1.2
}

finalize_execution() {
    local start_time=$1
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log "Total execution time: $((duration / 3600))h $(( (duration % 3600) / 60 ))m $((duration % 60))s"
    
    [[ -n "$EMAIL" ]] && {
        mail -s "Pipeline Completed: ${SPECIES_NAME}" "$EMAIL" < "$LOG_FILE"
        log "Notification email sent to: $EMAIL"
    }
}

# ============================
# Entry Point
# ============================
main "$@"
