#!/bin/bash
# UCE Probe Design Pipeline v2024.12.1
# Author: WJJ

# ========================
# Configuration
# ========================
set -euo pipefail  # Strict error handling
LOG_FILE="UCE_Probe_Design.log"
THREADS=8          # Adjust based on available cores

# ========================
# Initialize Logging
# ========================
{
    echo "=== Pipeline Started ==="
    echo "Timestamp: $(date)"
    echo "Parameters: $*"
} > "$LOG_FILE"

# ========================
# Runtime Tracking
# ========================
total_start=$(date +%s)

# ========================
# Pipeline Parameters
# ========================
declare -A STEP_DESCRIPTION=(
    [1]="FASTA to 2bit conversion"
    [2]="Reference genome setup"
    [3]="Read simulation"
    [4]="Read mapping"
    [5]="Probe table generation"
    [6]="Final probe design"
)

# ========================
# Function Definitions
# ========================

# Time-tracked execution wrapper
execute_step() {
    local step=$1
    local start end duration
    
    {
        echo -e "\n=== STEP $step: ${STEP_DESCRIPTION[$step]} ==="
        echo "Start: $(date)"
        start=$(date +%s)
    } | tee -a "$LOG_FILE"

    # Execute step commands
    case $step in
        1) process_fasta ;;
        2) setup_references ;;
        3) simulate_reads ;;
        4) map_reads ;;
        5) generate_tables ;;
        6) final_design ;;
    esac

    end=$(date +%s)
    duration=$((end - start))
    {
        echo "End: $(date)"
        echo "Duration: $duration seconds"
        echo "=== STEP $step COMPLETED ==="
    } | tee -a "$LOG_FILE"
}

process_fasta() {
    cd ./0-genomes || die "Missing 0-genomes directory"
    parallel -j $THREADS -0 "faToTwoBit {} {.}.2bit && mkdir -p {.} && mv {.}* {.}/" :::: <(find . -name "*.fasta" -print0)
}

setup_references() {
    mkdir -p {1-base_genome,2-inner_genomes,3-out_genome}
    cp ./0-genomes/{Basal*,Inner*,Outer*}* ./{1-base_genome,2-inner_genomes,3-out_genome}
    
    cd 1-base_genome || die "Missing base genome directory"
    stampy.py -G "$base" "$base.fasta" && stampy.py -g "$base" -H "$base"
}

simulate_reads() {
    cd ./2-inner_genomes || die "Missing inner genomes directory"
    parallel -j $THREADS "art_illumina --paired --in {} --out {.}-pe100-reads --len 100 --fcov 2" ::: *.fasta
}

map_reads() {
    mkdir -p 4-mapping && cd 4-mapping || die "Cannot create mapping directory"
    parallel -j 3 "stampy.py -g ../1-base_genome/$base -h ../1-base_genome/$base -M {} | samtools view -@ $THREADS -Sb - > {.}.bam" ::: ../2-inner_genomes/*-pe100-reads.fq.gz
}

generate_tables() {
    cd ./4-mapping || die "Missing mapping directory"
    phyluce_probe_get_multi_merge_table --conf bed_files.conf --base-taxon "$base" --output "${base}Temporary_share_count.bed"
    # Additional phyluce processing commands...
}

final_design() {
    mkdir -p 5-fina_probe_design && cd 5-fina_probe_design || die "Cannot create final design directory"
    # Final probe design commands...
}

# ========================
# Main Execution
# ========================
validate_input "$@"
setup_directories

for step in $(echo "${1:-123456}" | grep -o .); do
    [[ $step =~ [1-6] ]] || die "Invalid step: $step"
    execute_step "$step"
done

# ========================
# Final Reporting
# ========================
total_end=$(date +%s)
total_duration=$((total_end - total_start))

{
    echo -e "\n=== Pipeline Completed ==="
    echo "Total runtime: $total_duration seconds"
    echo "Breakdown:"
    for step in "${!STEP_DESCRIPTION[@]}"; do
        echo " - Step $step: ${STEP_DESCRIPTION[$step]} - ${DURATIONS[$step]:-0}s"
    done
} | tee -a "$LOG_FILE"
