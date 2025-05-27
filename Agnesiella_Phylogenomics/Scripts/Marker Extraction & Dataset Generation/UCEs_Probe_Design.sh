#!/bin/bash
# Script Name: UCE_Probe_Design_Workflow.sh
# Author: WJJ
# Date: 2024.12.30
# Description: This script automates the process of designing UCE probes using the PHYLUCE pipeline.
#              It involves genome preparation, read simulation, mapping, and probe design steps.

# --- Configuration Variables ---
# All user-configurable parameters are defined here.
# Please modify these variables according to your specific project.

# Prefix for the reference genome file (e.g., "reference_species")
# Example: BASE_GENOME_PREFIX="Homo_sapiens"
BASE_GENOME_PREFIX="your_base_genome_prefix"

# Space-separated list of prefixes for inner group genome files (e.g., "speciesA speciesB speciesC")
# Example: INNER_GENOMES_PREFIXES="Pan_troglodytes Gorilla_gorilla Pongo_abelii"
INNER_GENOMES_PREFIXES="your_inner_genome_prefix1 your_inner_genome_prefix2"

# Prefix for the outgroup genome file (e.g., "outgroup_species")
# Example: OUTGROUP_GENOME_PREFIX="Mus_musculus"
OUTGROUP_GENOME_PREFIX="your_outgroup_genome_prefix"

# Name of your taxonomic group/clade (e.g., "Primates")
GROUP_NAME="your_group_name"

# Initial minimum number of genomes that must share a UCE locus for it to be considered.
# This is used in Step 5 to select candidate UCEs.
INITIAL_SHARED_UCE_COUNT="initial_shared_uce_loci_count" # e.g., 4

# Final minimum number of genomes that must share a UCE locus for final probe design.
# This is used in Step 6 to refine the probe set.
FINAL_SHARED_UCE_COUNT="final_shared_uce_loci_count" # e.g., 8

# Number of CPU cores to use for parallel processing where applicable.
# Adjust based on your system's capabilities.
NUM_CORES_FATOTWOBIT=8
NUM_CORES_ART_ILLUMINA=8
NUM_CORES_STAMPY_MAP=3 # Note: STAMPY_MAP itself uses -t 12, but parallel jobs are limited to 3.
NUM_CORES_LASTZ_SQLITE=20 # For phyluce_probe_run_multiple_lastzs_sqlite

# --- Script Setup ---
# Exit immediately if a command exits with a non-zero status.
# Exit if an unset variable is used.
# Exit if any command in a pipeline fails.
set -euo pipefail

# Get the directory where the script is located and change to it.
# This ensures all relative paths work correctly regardless of where the script is called from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "$SCRIPT_DIR"

# --- Functions ---

# Function to split a number string (e.g., "123") into an array of individual digits ("1" "2" "3").
# This is used to control which steps of the workflow to execute.
# Usage: split_number "123" -> returns "1 2 3"
function split_number {
    local number_string="$1"
    local result="${number_string:0:1}"
    for (( i=1; i<${#number_string}; i++ )); do
        result+=" ${number_string:i:1}"
    done
    echo "$result"
}

# Convert FASTA genome file to 2bit format.
# Requires 'phyluce' conda environment.
# Args: $1 = path to genome FASTA file
FATOTWOBIT() {
    eval "$(conda shell.bash hook)"
    conda activate phyluce || { echo "Error: Could not activate phyluce conda environment. Please ensure it's installed."; exit 1; }
    
    local genome_file="$1"
    local genome_prefix=$(basename "$genome_file" ".fasta")
    
    echo "  Converting $genome_file to 2bit format..."
    faToTwoBit "$genome_file" "${genome_prefix}.2bit"
    
    # Create a subdirectory for the genome and move its files there
    mkdir -p "$genome_prefix"
    mv "${genome_prefix}"* "$genome_prefix"/
    echo "  Finished converting $genome_file."
}
export -f FATOTWOBIT # Export function for use with 'parallel'

# Simulate Illumina paired-end reads from a genome file using ART.
# Requires 'phyluce' conda environment.
# Args: $1 = path to genome FASTA file
ART_ILLUMINA() {
    eval "$(conda shell.bash hook)"
    conda activate phyluce || { echo "Error: Could not activate phyluce conda environment. Please ensure it's installed."; exit 1; }
    
    local genome_file="$1"
    local genome_prefix=$(basename "$genome_file" ".fasta")

    echo "  Simulating Illumina reads for $genome_file..."
    art_illumina --paired --in "${genome_file}" --out "${genome_prefix}-pe100-reads" \
                 --len 100 --fcov 2 --mflen 200 --sdev 150 \
                 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
    
    echo "  Reformatting simulated reads..."
    reformat.sh in1="${genome_prefix}-pe100-reads1.fq" in2="${genome_prefix}-pe100-reads2.fq" \
                out="${genome_prefix}-pe100-reads.fq.gz"
    
    # Clean up intermediate files
    rm "${genome_prefix}-pe100-reads1.fq" "${genome_prefix}-pe100-reads2.fq"
    echo "  Finished simulating reads for $genome_file."
}
export -f ART_ILLUMINA # Export function for use with 'parallel'

# Map simulated reads to the base genome using Stampy and process results.
# Requires 'stampy' and 'phyluce' conda environments.
# Args: $1 = path to gzipped FASTQ file of simulated reads
STAMPY_MAP() {
    eval "$(conda shell.bash hook)"
    conda activate stampy || { echo "Error: Could not activate stampy conda environment. Please ensure it's installed."; exit 1; }
    
    local genome_file="$1" # This is the simulated reads file, e.g., "speciesA-pe100-reads.fq.gz"
    local genome_prefix=$(basename "$genome_file" "-pe100-reads.fq.gz")
    
    echo "  Mapping reads from $genome_file to base genome ($BASE_GENOME_PREFIX)..."
    # Original script used $base for -g and -h, which is BASE_GENOME_PREFIX
    stampy.py --maxbasequal 93 -g "../1-base_genome/$BASE_GENOME_PREFIX" \
              -h "../1-base_genome/$BASE_GENOME_PREFIX" --substitutionrate=0.05 -t 12 \
              --insertsize=200 -M "$genome_file" | \
    samtools view -Sb -@ 12 - > "${genome_prefix}-to-${BASE_GENOME_PREFIX}.bam"
    
    echo "  Processing BAM file and generating BED file..."
    samtools view -h -F 4 -b -@ 12 "${genome_prefix}-to-${BASE_GENOME_PREFIX}.bam" | \
    bedtools bamtobed -i - -bed12 | \
    bedtools sort -i - | \
    bedtools merge -i - > "${genome_prefix}-to-${BASE_GENOME_PREFIX}-MAPPING.sort.merge.bed"
    
    # Activate phyluce for the next step, as it's used by phyluce_probe_strip_masked_loci_from_set
    conda activate phyluce || { echo "Error: Could not activate phyluce conda environment. Please ensure it's installed."; exit 1; }
    
    echo "  Stripping masked loci from BED file..."
    phyluce_probe_strip_masked_loci_from_set \
        --bed "${genome_prefix}-to-${BASE_GENOME_PREFIX}-MAPPING.sort.merge.bed" \
        --twobit "../1-base_genome/${BASE_GENOME_PREFIX}.2bit" \
        --output "${genome_prefix}-to-${BASE_GENOME_PREFIX}-MAPPING.sort.merge.strip.bed" \
        --filter-mask 0.25 --min-length 80
    
    # Append the path to the stripped BED file to the configuration file
    echo "${genome_prefix}:${genome_prefix}-to-${BASE_GENOME_PREFIX}-MAPPING.sort.merge.strip.bed" >> bed_files.conf
    echo "  Finished mapping and processing for $genome_file."
}
export -f STAMPY_MAP # Export function for use with 'parallel'

# --- Main Workflow ---

# Get the starting steps from the first command-line argument.
# Example: If you run `./script.sh 124`, it will execute steps 1, 2, and 4.
# If no argument is provided, it defaults to executing all steps (1 through 6).
STARTING_STEP=($(split_number "${1:-123456}"))

# Create species.list if it doesn't exist.
# This file lists all genome prefixes found in the 0-genomes directory.
if [ ! -f "species.list" ]; then
    echo "Creating species.list from 0-genomes directory..."
    # Ensure 0-genomes directory exists before trying to list files
    if [ ! -d "./0-genomes" ]; then
        echo "Error: Directory ./0-genomes not found. Please place your FASTA genomes here."
        exit 1
    fi
    find "./0-genomes" -maxdepth 1 -name "*.fasta" -type f -print0 | \
        xargs -0 -n 1 basename -s ".fasta" > species.list
    echo "species.list created."
fi

# Loop through the specified steps
for num in "${STARTING_STEP[@]}"
do
    echo "" # Newline for readability
    echo "--- Starting Step $num ---"
    case $num in
        1)
            echo "Step 1: Converting all genomes in 0-genomes to 2bit format."
            # Run faToTwoBit in parallel for all FASTA files in the 0-genomes directory
            find "./0-genomes" -name "*.fasta" -type f -print0 | \
                parallel -0 -n 1 -I {} -j "$NUM_CORES_FATOTWOBIT" FATOTWOBIT {}
            echo "Step 1 completed."
            ;;
        2)
            echo "Step 2: Setting up genome directories and indexing the base genome."
            # Create necessary directories
            mkdir -p "1-base_genome" "2-inner_genomes" "3-out_genome"
            
            # Copy genomes (FASTA and 2bit files) to their respective directories
            echo "  Copying base genome files to 1-base_genome..."
            cp "./0-genomes/${BASE_GENOME_PREFIX}/${BASE_GENOME_PREFIX}.fasta" "./1-base_genome/"
            cp "./0-genomes/${BASE_GENOME_PREFIX}/${BASE_GENOME_PREFIX}.2bit" "./1-base_genome/"
            
            echo "  Copying inner genomes files to 2-inner_genomes..."
            for i in $INNER_GENOMES_PREFIXES; do
                cp "./0-genomes/${i}/${i}.fasta" "./2-inner_genomes/"
                cp "./0-genomes/${i}/${i}.2bit" "./2-inner_genomes/"
            done
            
            echo "  Copying outgroup genome files to 3-out_genome..."
            cp "./0-genomes/${OUTGROUP_GENOME_PREFIX}/${OUTGROUP_GENOME_PREFIX}.fasta" "./3-out_genome/"
            cp "./0-genomes/${OUTGROUP_GENOME_PREFIX}/${OUTGROUP_GENOME_PREFIX}.2bit" "./3-out_genome/"
            
            # Activate stampy environment for indexing
            eval "$(conda shell.bash hook)"
            conda activate stampy || { echo "Error: Could not activate stampy conda environment. Please ensure it's installed."; exit 1; }
            
            echo "  Indexing base genome ($BASE_GENOME_PREFIX) for Stampy..."
            # Run Stampy indexing in a subshell to avoid changing the main script's directory
            (
                cd "1-base_genome" || { echo "Error: Could not change directory to 1-base_genome."; exit 1; }
                stampy.py --species="$BASE_GENOME_PREFIX" --assembly="$BASE_GENOME_PREFIX" -G "$BASE_GENOME_PREFIX" "${BASE_GENOME_PREFIX}.fasta"
                stampy.py -g "$BASE_GENOME_PREFIX" -H "$BASE_GENOME_PREFIX"
            )
            echo "Step 2 completed."
            ;;
        3)
            echo "Step 3: Simulating Illumina reads for inner genomes."
            # Ensure 2-inner_genomes directory exists and contains FASTA files
            if [ ! -d "./2-inner_genomes" ] || [ -z "$(find "./2-inner_genomes" -name "*.fasta" -print -quit)" ]; then
                echo "Error: Directory ./2-inner_genomes not found or empty. Please run step 2 first."
                exit 1
            fi

            # Run ART_ILLUMINA in parallel for all FASTA files in 2-inner_genomes
            find "./2-inner_genomes" -name "*.fasta" -type f -print0 | \
                parallel -0 -n 1 -I {} -j "$NUM_CORES_ART_ILLUMINA" ART_ILLUMINA {}
            echo "Step 3 completed."
            ;;
        4)
            echo "Step 4: Mapping simulated reads to the base genome and processing results."
            mkdir -p "4-mapping"
            
            # Clear or create a fresh bed_files.conf for this step.
            # This file will be populated by the STAMPY_MAP function with paths to stripped BED files.
            echo '[beds]' > "./4-mapping/bed_files.conf"
            
            # Export BASE_GENOME_PREFIX so STAMPY_MAP function can access it in subshells
            export BASE_GENOME_PREFIX
            
            # Run STAMPY_MAP in parallel for all simulated read files
            find "./2-inner_genomes" -name "*-pe100-reads.fq.gz" -type f -print0 | \
                parallel -0 -n 1 -I {} -j "$NUM_CORES_STAMPY_MAP" STAMPY_MAP {}
            echo "Step 4 completed."
            ;;
        5)
            echo "Step 5: Designing initial UCE probes based on shared loci."
            # Activate phyluce environment
            eval "$(conda shell.bash hook)"
            conda activate phyluce || { echo "Error: Could not activate phyluce conda environment. Please ensure it's installed."; exit 1; }
            
            # Run phyluce commands in a subshell from the 4-mapping directory
            (
                cd "4-mapping" || { echo "Error: Could not change directory to 4-mapping."; exit 1; }
                
                echo "  Getting multi-merge table from bed files..."
                phyluce_probe_get_multi_merge_table --conf bed_files.conf \
                                                    --base-taxon "$BASE_GENOME_PREFIX" \
                                                    --output "${BASE_GENOME_PREFIX}.sqlite"
                
                echo "  Querying multi-merge table for shared UCEs (initial count: $INITIAL_SHARED_UCE_COUNT)..."
                phyluce_probe_query_multi_merge_table --db "${BASE_GENOME_PREFIX}.sqlite" \
                                                      --base-taxon "$BASE_GENOME_PREFIX" \
                                                      --output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.bed" \
                                                      --specific-counts "$INITIAL_SHARED_UCE_COUNT"
                
                echo "  Getting genome sequences from BED file..."
                phyluce_probe_get_genome_sequences_from_bed \
                    --bed "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.bed" \
                    --twobit "../1-base_genome/${BASE_GENOME_PREFIX}.2bit" \
                    --buffer-to 160 \
                    --output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.fasta"
                
                echo "  Designing tiled probes..."
                phyluce_probe_get_tiled_probes \
                    --input "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.fasta" \
                    --probe-prefix "uce-" \
                    --design "${GROUP_NAME}-v1" \
                    --designer "wjj" \
                    --tiling-density 3 \
                    --overlap middle \
                    --masking 0.25 \
                    --remove-gc \
                    --two-probes \
                    --output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes"
                
                echo "  Running LASTZ for self-comparison of probes..."
                phyluce_probe_easy_lastz \
                    --target "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes" \
                    --query "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes" \
                    --identity 50 \
                    --coverage 50 \
                    --output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes-TO-SELF-PROBES.lastz"
                
                echo "  Removing duplicate hits from probes..."
                # Original command, no --output parameter added as per user request
                phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
                    --fasta "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes" \
                    --lastz "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp.probes-TO-SELF-PROBES.lastz" \
                    --probe-prefix="uce-"
            )
            echo "Step 5 completed."
            ;;
        6)
            echo "Step 6: Finalizing UCE probe design across all species."
            # Activate phyluce environment
            eval "$(conda shell.bash hook)"
            conda activate phyluce || { echo "Error: Could not activate phyluce conda environment. Please ensure it's installed."; exit 1; }
            
            mkdir -p "5-fina_probe_design"
            # Run phyluce commands in a subshell from the 5-fina_probe_design directory
            (
                cd "5-fina_probe_design" || { echo "Error: Could not change directory to 5-fina_probe_design."; exit 1; }
                
                # Create genome configuration file for phyluce_probe_run_multiple_lastzs_sqlite
                echo '[scaffolds]' > "${GROUP_NAME}-genome.conf"
                while IFS= read -r SPECIES; do
                    echo "$SPECIES:../0-genomes/$SPECIES/$SPECIES.2bit" >> "${GROUP_NAME}-genome.conf"
                done < "../species.list"
                
                echo "  Running multiple LASTZ alignments against all genomes..."
                phyluce_probe_run_multiple_lastzs_sqlite \
                    --probefile "../4-mapping/${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp-DUPE-SCREENED.probes" \
                    --scaffoldlist "$(cat ../species.list)" \
                    --genome-base-path "../0-genomes" \
                    --identity 50 \
                    --cores "$NUM_CORES_LASTZ_SQLITE" \
                    --db "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.sqlite" \
                    --output "${GROUP_NAME}-genome-lastz"
                
                echo "  Slicing sequences from genomes based on LASTZ results..."
                phyluce_probe_slice_sequence_from_genomes \
                    --conf "${GROUP_NAME}-genome.conf" \
                    --lastz "${GROUP_NAME}-genome-lastz" \
                    --probes 180 \
                    --name-pattern "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
                    --output "${GROUP_NAME}-genome-fasta"
                
                echo "  Getting multi-fasta table..."
                phyluce_probe_get_multi_fasta_table \
                    --fastas "${GROUP_NAME}-genome-fasta" \
                    --output "multifastas.sqlite" \
                    --base-taxon "$BASE_GENOME_PREFIX"
                
                echo "  Querying multi-fasta table for final shared UCEs (final count: $FINAL_SHARED_UCE_COUNT)..."
                phyluce_probe_query_multi_fasta_table \
                    --db "multifastas.sqlite" \
                    --base-taxon "$BASE_GENOME_PREFIX" \
                    --output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}-back-to-${FINAL_SHARED_UCE_COUNT}.conf" \
                    --specific-counts "$FINAL_SHARED_UCE_COUNT"
                
                echo "  Getting tiled probes from multiple inputs for final probe list..."
                # Original command, using --multi-fasta-output as per user request
                phyluce_probe_get_tiled_probe_from_multiple_inputs \
                    --fastas "${GROUP_NAME}-genome-fasta" \
                    --multi-fasta-output "${BASE_GENOME_PREFIX}+${INITIAL_SHARED_UCE_COUNT}-back-to-${FINAL_SHARED_UCE_COUNT}.conf" \
                    --probe-prefix "uce-" \
                    --designer "wjj" \
                    --design "${GROUP_NAME}-v1" \
                    --tiling-density 3 \
                    --overlap middle \
                    --masking 0.25 \
                    --remove-gc \
                    --two-probes \
                    --output "${GROUP_NAME}-v1-master-probe-list.fasta"
                
                echo "  Running LASTZ for self-comparison of final probes..."
                phyluce_probe_easy_lastz \
                    --target "${GROUP_NAME}-v1-master-probe-list.fasta" \
                    --query "${GROUP_NAME}-v1-master-probe-list.fasta" \
                    --identity 50 \
                    --coverage 50 \
                    --output "${GROUP_NAME}-v1-master-probe-list-TO-SELF-PROBES.lastz"
                
                echo "  Removing duplicate hits from final probes..."
                # Original command, no --output parameter added as per user request
                phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
                    --fasta "${GROUP_NAME}-v1-master-probe-list.fasta" \
                    --lastz "${GROUP_NAME}-v1-master-probe-list-TO-SELF-PROBES.lastz" \
                    --probe-prefix="uce-"
            )
            echo "Step 6 completed."
            ;;
        *)
            echo "Invalid step number: ${num}. Please enter a number within the range [1-6]."
            exit 1
            ;;
    esac
    echo "--- Finished Step $num ---"
done

echo "" # Newline for readability
echo "Workflow completed successfully!"
