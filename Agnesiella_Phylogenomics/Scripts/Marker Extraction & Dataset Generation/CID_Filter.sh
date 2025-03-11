#!/bin/bash
# CID Analysis Pipeline v2024.12.1
# Author: WJJ

# Usage:
#   bash CID_Filter.sh <loci_dir> <treefile_dir> <outgroup> <output_dir>
#
# Requirements:
#   - realpath, Rscript, python
#   - Calculate_CID_Distance.R and Visualize_CID_Results.py in PATH

set -euo pipefail

usage() {
    echo "Usage: $0 <loci_dir> <treefile_dir> <outgroup> <output_dir>"
    echo "  loci_dir:     Directory containing locus files"
    echo "  treefile_dir: Directory containing treefiles"
    echo "  outgroup:     Outgroup specification for analysis"
    echo "  output_dir:   Target directory for results"
    exit 1
}

validate_directory() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory '$1' not found"
        exit 1
    fi
}

# Argument validation
[ $# -eq 4 ] || usage

LOCI_DIR="$1"
TREEFILE_DIR="$2"
OUTGROUP="$3"
OUTPUT_DIR="$4"

validate_directory "$LOCI_DIR"
validate_directory "$TREEFILE_DIR"

# Create workspace
mkdir -p "$OUTPUT_DIR"
WORKSPACE=$(realpath "$OUTPUT_DIR")

# Setup directory structure
CID_DIRS=(
    "1-CID_loci"
    "2-CID_trees" 
    "3-CID_remove_loci"
    "4-CID_remove_trees"
)

for dir in "${CID_DIRS[@]}"; do
    mkdir -p "${WORKSPACE}/${dir}"
done

# Copy input data
cp -r "${LOCI_DIR}/"* "${WORKSPACE}/1-CID_loci/"
cp -r "${TREEFILE_DIR}/"* "${WORKSPACE}/2-CID_trees/"

cd "$WORKSPACE"

# Generate file lists
TREEFILE_COUNT=$(find 2-CID_trees -name '*.treefile' | wc -l)
[ $TREEFILE_COUNT -eq 0 ] && { echo "Error: No treefiles found"; exit 1; }

cat 2-CID_trees/*.treefile > 1-ABS_genetrees.treefile
find 1-CID_loci -maxdepth 1 -type f > 2-ABS_genetrees.list

# Dependency checks
command -v Rscript >/dev/null || { echo "Error: Rscript required"; exit 1; }
command -v python >/dev/null || { echo "Error: Python required"; exit 1; }

# Execute analysis
Rscript Calculate_CID_Distance.R "$WORKSPACE" 1-ABS_genetrees.treefile || {
    echo "R script failed"; exit 1;
}

python Visualize_CID_Results.py 3-CID_distance_matrix.csv 2-ABS_genetrees.list || {
    echo "Python script failed"; exit 1;
}

# Handle outliers
REMOVE_LIST="5-Exceed_MAD_CID_trees.txt"
[ -f "$REMOVE_LIST" ] && {
    while read -r target; do
        mv -f "1-CID_loci/${target}" "3-CID_remove_loci/" || true
        mv -f "2-CID_trees/${target}.treefile" "4-CID_remove_trees/" || true
    done < "$REMOVE_LIST"
}

echo "Processing completed. Results in: ${WORKSPACE}"
