#!/bin/bash
# UCE Extraction Pipeline v2024.12.1
# Author: WJJ

set -euo pipefail

# ========================
# Global Variables
# ========================
start_time=$(date +%s)
declare -A step_times=()
THREADS=8
LOG_FILE="uce_pipeline.log"

# ========================
# Function Definitions
# ========================

log() {
    local message="$1"
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $message" | tee -a "$LOG_FILE"
}

record_step_time() {
    local step_name="$1"
    local step_start="$2"
    local step_end=$(date +%s)
    step_times["$step_name"]=$((step_end - step_start))
}

validate_inputs() {
    local start=$(date +%s)
    
    if [ "$#" -ne 2 ]; then
        log "Error: Requires 2 arguments"
        echo "Usage: $0 <probe_fasta> <assembly_directory>"
        exit 1
    fi

    INPUT_PROBE=$(realpath "${1//\'/}")
    DIR_ASSEMBLY=$(realpath "${2//\'/}")

    [ -f "$INPUT_PROBE" ] || { log "Probe file missing"; exit 1; }
    [ -d "$DIR_ASSEMBLY" ] || { log "Assembly dir missing"; exit 1; }

    record_step_time "Input Validation" "$start"
}

setup_directories() {
    local start=$(date +%s)
    
    log "Creating directory structure..."
    mkdir -p {0-genomes,1-uce_align,2-uces,3-raw_loci,4-loci_filter}
    
    record_step_time "Directory Setup" "$start"
}

# ========================
# Main Pipeline
# ========================
{
    log "=== Starting UCE Pipeline ==="
    
    # Phase 1: Initialization
    validate_inputs "$@"
    setup_directories

    # Phase 2: Data Preparation
    prep_start=$(date +%s)
    log "Copying probe set..."
    cp "$INPUT_PROBE" probe.fasta
    log "Transferring assemblies..."
    cp -r "$DIR_ASSEMBLY"/* 0-genomes/
    record_step_time "Data Preparation" "$prep_start"

    # Phase 3: Genome Processing
    process_start=$(date +%s)
    (
        cd 0-genomes/
        log "Standardizing genome names..."
        for dir in */; do
            base=$(basename "$dir")
            first_char=${base:0:1}
            rest=${base:1}
            [[ "$first_char" =~ [a-z] ]] && mv "$base" "$(echo "$first_char" | tr 'a-z' 'A-Z')$rest"
        done
    )
    record_step_time "Genome Processing" "$process_start"

    # Phase 4: Alignment Operations
    align_start=$(date +%s)
    log "Running LASTZ alignments..."
    phyluce_probe_run_multiple_lastzs_sqlite \
        --db 1-uce_align/uces.sqlite \
        --output 1-uce_align/target-genome-lastz \
        --probefile probe.fasta \
        --scaffoldlist species.list \
        --genome-base-path . \
        --identity 50 \
        --cores "$THREADS"
    record_step_time "Sequence Alignment" "$align_start"

    # Additional processing steps would follow similar patterns...

    # Final Reporting
    end_time=$(date +%s)
    total_time=$((end_time - start_time))
    
    log "\n=== Time Summary ==="
    for step in "${!step_times[@]}"; do
        log "${step}: ${step_times[$step]} seconds"
    done
    log "Total execution time: $total_time seconds"

} 2>&1 | tee -a "$LOG_FILE"
