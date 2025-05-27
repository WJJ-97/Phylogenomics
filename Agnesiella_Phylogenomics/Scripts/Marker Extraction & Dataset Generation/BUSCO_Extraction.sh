#!/bin/bash
#2021.08.04 by ZF
#2023.06.28 by DSY
#2024.12.31 by WJJ

# Script: BUSCO Single-Copy Loci Extractor and Filter
# Description: This script processes BUSCO results to extract single-copy orthologs,
#              merges them by locus, and filters out loci with insufficient taxa.

set -euo pipefail

# --- Usage ---
# This script processes BUSCO results from multiple species to extract,
# merge, and filter single-copy BUSCO sequences.
#
# Usage: ./your_script_name.sh
#
# It will prompt for input and output directories.
# Input directory should contain subdirectories for each species,
# where each species directory contains the BUSCO run results (e.g., run_lineage).
# Output directory will be created to store processed files.

echo "--- BUSCO Sequence Processing Script ---"

# input the name of input directory
read -p "Please input the name of input directory containing all BUSCO results folders (e.g. /home/zf/Desktop/materials/datasets_examples/BUSCOs): " DIR_INPUT_TEMP
# Use realpath for absolute path and remove potential quotes
DIR_INPUT=$(realpath "$(echo "$DIR_INPUT_TEMP" | sed "s/'//g")")

# Check if input directory exists
if [ ! -d "$DIR_INPUT" ]; then
    echo "Error: Input directory '$DIR_INPUT' not found." >&2
    exit 1
fi
echo "Input directory set to: $DIR_INPUT"

# input the name of output directory
read -p "Please input the name of output directory, or an existing directory: " DIR_OUTPUT_TEMP
# Construct output directory path relative to PWD and remove potential quotes
DIR_OUTPUT="${PWD}/$(echo "$DIR_OUTPUT_TEMP" | sed "s/'//g")"

# Create output directory (use -p to avoid error if it exists)
mkdir -p "$DIR_OUTPUT"
echo "Output directory set to: $DIR_OUTPUT"

# --- Phase 1: Collect and Prepare Single-Copy BUSCO Sequences ---
echo "--- Phase 1: Collecting and preparing sequences ---"
DIR_SINGLE_COPY="$DIR_OUTPUT/1-single_copy_busco_sequences"
mkdir -p "$DIR_SINGLE_COPY"

# Generate a file containing species list using find for robustness
# Find directories one level deep in the input directory
find "$DIR_INPUT" -maxdepth 1 -minddepth 1 -type d -printf '%f\n' | sort > "$DIR_OUTPUT/species.list"

# Check if species list is empty
if [ ! -s "$DIR_OUTPUT/species.list" ]; then
    echo "Error: No species directories found in '$DIR_INPUT'." >&2
    exit 1
fi

echo "Species list generated: $(cat "$DIR_OUTPUT/species.list" | wc -l) species found."

# Temporary file to collect all unique locus names
TEMP_ALL_LOCI="$DIR_SINGLE_COPY/temp_all_loci.list"
> "$TEMP_ALL_LOCI" # Create/clear the temp file

# Read species list line by line
while IFS= read -r SPECIES; do
    echo "Processing species: $SPECIES"
    SPECIES_SEQ_DIR="$DIR_SINGLE_COPY/${SPECIES}_sequence"
    mkdir -p "$SPECIES_SEQ_DIR"

    # Find BUSCO run directory (e.g., run_lineage)
    BUSCO_RUN_DIR=$(find "$DIR_INPUT/$SPECIES" -maxdepth 1 -type d -name "run_*" | head -n 1)

    if [ -z "$BUSCO_RUN_DIR" ]; then
        echo "Warning: No BUSCO run directory found for species '$SPECIES'. Skipping." >&2
        continue # Skip to the next species
    fi

    BUSCO_SEQ_SOURCE="$BUSCO_RUN_DIR/busco_sequences/single_copy_busco_sequences"

    if [ ! -d "$BUSCO_SEQ_SOURCE" ]; then
         echo "Warning: Single copy sequences directory not found for species '$SPECIES'. Skipping." >&2
         continue # Skip to the next species
    fi

    # Copy and modify the head name of the fasta files for each locus
    # Use find with -print0 and while read -d $'\0' for robust file handling
    find "$BUSCO_SEQ_SOURCE" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.faa" \) -print0 | while IFS= read -r -d $'\0' seq_file; do
        # Get filename without path
        filename=$(basename "$seq_file")
        # Copy the file
        cp "$seq_file" "$SPECIES_SEQ_DIR/$filename"

        # Modify the header line (starts with >) to just the species name
        # Use sed -i for in-place editing
        sed -i "s/^>.*/>${SPECIES}/" "$SPECIES_SEQ_DIR/$filename"

        # --- Retained as requested ---
        # Remove all asterisks from the file. Note: In FAA files, * represents stop codons.
        # Removing them might affect downstream protein sequence analysis.
        sed -i "s/*//g" "$SPECIES_SEQ_DIR/$filename"
        # --- End Retained ---

        # Extract locus name (filename without extension) and add to temp list
        locus_name="${filename%.*}" # Remove extension (.fna or .faa)
        echo "$locus_name" >> "$TEMP_ALL_LOCI"

    done

done < "$DIR_OUTPUT/species.list" # Read species list into the loop

# Generate a file of unique loci list from the temporary list
sort "$TEMP_ALL_LOCI" | uniq > "$DIR_SINGLE_COPY/loci.list"
rm -f "$TEMP_ALL_LOCI" # Clean up temp file

echo "Unique loci list generated: $(cat "$DIR_SINGLE_COPY/loci.list" | wc -l) loci found."

# --- Phase 2: Merge sequences of the same locus into raw files ---
echo "--- Phase 2: Merging sequences by locus ---"
DIR_RAW_LOCI="$DIR_OUTPUT/2-raw_loci"
mkdir -p "$DIR_RAW_LOCI/fna" "$DIR_RAW_LOCI/faa"
ABSENCE_LOG="$DIR_RAW_LOCI/absence.log"
> "$ABSENCE_LOG" # Create/clear absence log

# Read loci list line by line
while IFS= read -r LOCI; do
    # Read species list line by line for each locus
    while IFS= read -r SPECIES; do
        FNA_FILE="$DIR_SINGLE_COPY/${SPECIES}_sequence/${LOCI}.fna"
        FAA_FILE="$DIR_SINGLE_COPY/${SPECIES}_sequence/${LOCI}.faa"

        # Check if the locus file exists for this species
        if [ -f "$FNA_FILE" ]; then
            # Append sequences to the merged locus file
            cat "$FNA_FILE" >> "$DIR_RAW_LOCI/fna/$LOCI.fna"
            cat "$FAA_FILE" >> "$DIR_RAW_LOCI/faa/$LOCI.faa"
        else
            # Log absence
            echo "$LOCI in $SPECIES does not exist" >> "$ABSENCE_LOG"
        fi
    done < "$DIR_OUTPUT/species.list" # Read species list into the inner loop
done < "$DIR_SINGLE_COPY/loci.list" # Read loci list into the outer loop

echo "Sequence merging complete. Absence log generated: $ABSENCE_LOG"

# --- Phase 3: Filter loci having too few taxa (less than three) ---
echo "--- Phase 3: Filtering loci ---"
DIR_FILTERED_LOCI="$DIR_OUTPUT/3-loci_filter"
mkdir -p "$DIR_FILTERED_LOCI/fna" "$DIR_FILTERED_LOCI/faa"
LOCI_REMAIN_LIST="$DIR_FILTERED_LOCI/loci_name_remain.list"
> "$LOCI_REMAIN_LIST" # Create/clear remaining loci list

MIN_PRESENT_TAXA=3

echo "Filtering loci: Keeping loci present in at least $MIN_PRESENT_TAXA species."

# Iterate through the merged FNA files in the raw directory
# Use find with -print0 for robust file handling
find "$DIR_RAW_LOCI/fna" -maxdepth 1 -type f -name "*.fna" -print0 | while IFS= read -r -d $'\0' raw_fna_file; do
    # Get locus name from filename
    LOCI=$(basename "$raw_fna_file" .fna)
    raw_faa_file="$DIR_RAW_LOCI/faa/$LOCI.faa"

    # Count the number of sequences (headers starting with >) in the merged file
    HEADER_COUNT=$(grep -c "^>" "$raw_fna_file")

    # Check if the number of sequences meets the minimum requirement
    if [ "$HEADER_COUNT" -ge "$MIN_PRESENT_TAXA" ]; then
        # Copy the locus files to the filtered directory
        cp "$raw_fna_file" "$DIR_FILTERED_LOCI/fna/"
        cp "$raw_faa_file" "$DIR_FILTERED_LOCI/faa/"
        # Log the name of the remaining locus
        echo "$LOCI" >> "$LOCI_REMAIN_LIST"
    fi
done

echo "Loci filtering complete. Remaining loci list: $LOCI_REMAIN_LIST"
echo "Number of remaining loci: $(cat "$LOCI_REMAIN_LIST" | wc -l)"

# --- Phase 4: Rename files and count ---
echo "--- Phase 4: Renaming files ---"

# Rename .fna to .fas in the filtered directory using bash loop
find "$DIR_FILTERED_LOCI/fna" -maxdepth 1 -type f -name "*.fna" -print0 | while IFS= read -r -d $'\0' fna_file; do
    mv "$fna_file" "${fna_file%.fna}.fas"
done

# Rename .faa to .fas in the filtered directory using bash loop
find "$DIR_FILTERED_LOCI/faa" -maxdepth 1 -type f -name "*.faa" -print0 | while IFS= read -r -d $'\0' faa_file; do
    mv "$faa_file" "${faa_file%.faa}.fas"
done

echo "Renaming complete (.fna -> .fas, .faa -> .fas)."

# Count the number of remaining FAA/FAS files
FINAL_COUNT=$(find "$DIR_FILTERED_LOCI/faa" -maxdepth 1 -type f -name "*.fas" | wc -l)
echo "Final count of filtered FAA/FAS files: $FINAL_COUNT"

echo "--- Script finished successfully ---"
