#!/usr/bin/env bash
# BUSCO Extraction v2024.12.1
# Author: WJJ

# Initialize performance tracking
SECONDS=0
declare -A timing=()

# Configuration
THREADS_PER_PHASE=2  # Threads allocation per processing stage

{{
# Helper function for time logging
log_time() {
    local stage_name=$1
    timing["$stage_name"]=$SECONDS
    echo "===== [TIMING] $stage_name: $(date +'%T')" | tee -a "$LOG_FILE"
    SECONDS=0
}

# Main processing stages
initialize_processing() {
    echo "===== INITIALIZATION STARTED ====="
    # Input validation and setup
    [[ $# -lt 2 ]] && { echo "Error: Insufficient arguments"; usage; exit 1; }
    validate_dependencies
    setup_directories
    log_time "Initialization"
}

process_species_data() {
    echo "===== SPECIES PROCESSING STARTED ====="
    generate_species_list
    parallel_extract_sequences
    compile_loci_list
    log_time "SpeciesProcessing"
}

merge_loci_data() {
    echo "===== LOCI MERGING STARTED ====="
    parallel_merge_loci
    log_time "LociMerging"
}

filter_loci_data() {
    echo "===== LOCI FILTERING STARTED ====="
    calculate_thresholds
    parallel_filter_loci
    finalize_output
    log_time "LociFiltering"
}

# Enhanced dependency check
validate_dependencies() {
    local deps=("parallel" "awk" "sed")
    for cmd in "${deps[@]}"; do
        command -v "$cmd" >/dev/null 2>&1 || {
            echo "Error: Required command '$cmd' not found"
            exit 1
        }
    done
}

# Optimized directory setup
setup_directories() {
    INPUT_DIR=$(realpath "${1//\'/}")
    OUTPUT_DIR=$(realpath -m "${2//\'/}")
    LOG_FILE="$OUTPUT_DIR/processing.log"
    
    mkdir -p "$OUTPUT_DIR"/{1-single_copy,2-raw_loci,3-filtered_loci}
    exec 3>&1 4>&2
    exec > >(tee -a "$LOG_FILE") 2>&1
}

# Main execution flow
main() {
    initialize_processing "$@"
    
    process_species_data
    merge_loci_data
    filter_loci_data
    
    # Performance summary
    echo "===== PROCESSING COMPLETED ====="
    printf "Stage breakdown:\n"
    printf "%-15s %s\n" "Stage" "Seconds"
    for stage in "${!timing[@]}"; do
        printf "%-15s %d\n" "$stage" "${timing[$stage]}"
    done
    printf "Total runtime: %d seconds\n" $SECONDS
}

# Execution entry
main "$@"
}}

# Implementation details for key functions:

generate_species_list() {
    find "$INPUT_DIR" -maxdepth 1 -type d -printf '%f\n' | grep -v "^$" > "$OUTPUT_DIR/species.list"
}

parallel_extract_sequences() {
    local species_count=$(wc -l < "$OUTPUT_DIR/species.list")
    local thread_count=$(( species_count < THREADS_PER_PHASE ? species_count : THREADS_PER_PHASE ))
    
    parallel -j "$thread_count" --bar --eta '
    species={}
    seq_dir="$OUTPUT_DIR/1-single_copy/${species}_seq"
    mkdir -p "$seq_dir"
    
    # Fast header modification using AWK
    find "$INPUT_DIR/$species" -name "*.f?a" | parallel -j 2 '
        awk -v sp="$species" '"'"'
        /^>/ {print ">"sp; next} 
        {print}
        '"'"' {} > "$seq_dir/$(basename {})"
    '
    ' :::: "$OUTPUT_DIR/species.list"
}

compile_loci_list() {
    find "$OUTPUT_DIR/1-single_copy" -name "*.f?a" -printf '%f\n' | 
    cut -d. -f1 | sort -u > "$OUTPUT_DIR/loci.masterlist"
}
