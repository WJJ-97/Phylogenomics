#!/bin/bash
#2020.08.04 by ZF
#2023.07.11 by DSY
#2024.12.30 by WJJ

# De Novo Genome Assembly and Assessment Script
#
# Required Tools: pigz, fastp, BBtools (bbnorm.sh), Redundans (redundans.py),
#                 Minimap2, Samtools, BESST (runBESST), GapCloser, SPAdes (spades.py),
#                 SeqKit, BUSCO
#
# Usage:
#   bash <script_name> <read1.fq.gz> <read2.fq.gz> <species_name> <steps_to_run>
#
# Example:
#   bash assemble_genome.sh reads_R1.fq.gz reads_R2.fq.gz MySpecies 12345678
#

# Exit immediately if a command exits with a non-zero status.
set -e

# --- 0. Initialization and Argument Parsing ---

start_time=$(date +%s) # Record script start time

# Check number of command-line arguments
# Now expecting 4 arguments instead of 5
if [ "$#" -ne 4 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <species_name> <steps_to_run>"
    echo "Example: $0 reads_R1.fq.gz reads_R2.fq.gz MySpecies 12345678"
    exit 1
fi

# Assign command-line arguments to descriptive variables
INPUT_R1="$1"
INPUT_R2="$2"
SPECIES="$3"
STARTING_STEPS_RAW="$4" # e.g., "12345678"

# Extract base names and prefixes from input filenames
filename1=$(basename "${INPUT_R1}")
filename2=$(basename "${INPUT_R2}")
prefix1="${filename1%%.*}" # e.g., "reads_R1" from "reads_R1.fq.gz"
prefix2="${filename2%%.*}" # e.g., "reads_R2" from "reads_R2.fq.gz"

# --- 1. Define Global Parameters ---
READ_LENGTH=150 # Read length of sequencing data
THREADS=16      # Number of threads to use
MEMORY_MB=200000 # Total memory limit (MB), e.g., 200GB
BBTOOLS_MEM_GB="80g" # Specific memory for BBtools (GB)
SAMTOOLS_MEM_PER_THREAD="1706M" # Samtools sort memory per thread (MB)

# --- IMPORTANT: BUSCO Lineage Path (Hardcoded) ---
# Please set the absolute path to your BUSCO lineage dataset here.
# Example: /path/to/busco_lineage/eukaryota_odb10
DIR_BUSCO_LINEAGE="/path/to/your/busco_lineage/eukaryota_odb10" # <--- SET THIS PATH

# --- 2. Set Up Working Directories and Log File ---
LOG_DIR="0-log"
BUSCO_SUM_DIR="1-busco_sum"
SPECIES_BASE_DIR="${SPECIES}" # Main output directory for species-specific results

# Create necessary directories, -p ensures parent directories are also created
mkdir -p "${LOG_DIR}"
mkdir -p "${BUSCO_SUM_DIR}"
mkdir -p "${SPECIES_BASE_DIR}"

# Define the absolute path for the main log file
LOG_FILE="$(realpath "${LOG_DIR}/${SPECIES}_${STARTING_STEPS_RAW}_assemble.log")"

# Redirect all subsequent stdout and stderr to the log file, and also display on terminal
# Note: After this line, all echo and command output will automatically be written to LOG_FILE
exec &> >(tee -a "${LOG_FILE}")

echo "--- Genome Assembly Script Started ---"
echo "Date: $(date)"
echo "Input File R1: ${INPUT_R1}"
echo "Input File R2: ${INPUT_R2}"
echo "Species Name: ${SPECIES}"
echo "Steps to Run: ${STARTING_STEPS_RAW}"
echo "BUSCO Lineage Path (Hardcoded): ${DIR_BUSCO_LINEAGE}" # Indicate it's hardcoded
echo "Read Length: ${READ_LENGTH}"
echo "Threads: ${THREADS}"
echo "Memory (MB): ${MEMORY_MB}"
echo "--------------------------------------"
echo ""

# --- 3. Dependency Tool Check Function ---
# Function: Checks if a tool exists and sets its executable path
# Arguments:
#   $1: Tool name (e.g., "pigz")
#   $2: Variable name to store executable path (e.g., "PIGZ_EXE")
#   $3: Example path hint for user input (e.g., "/usr/bin")
check_tool() {
    local tool_name="$1"
    local var_name="$2"
    local example_path_hint="$3"
    local exe_path=""
    local dir_path=""

    echo -n "Checking ${tool_name} ...... "

    # First, try to find the tool in the PATH environment variable
    exe_path=$(which "${tool_name}" 2>/dev/null)

    if [ -x "${exe_path}" ]; then
        echo "OK"
        # Set global variable, e.g., PIGZ_EXE="/usr/bin/pigz"
        eval "${var_name}='${exe_path}'"
    else
        echo "Not Found"
        # Loop to prompt user for installation directory until executable is found
        until [ -x "${dir_path}/${tool_name}" ]; do
            read -p "${tool_name} not found. Please input its absolute installation directory (e.g., ${example_path_hint}): " dir_path
            # Remove trailing slash from path (if any)
            dir_path="${dir_path%/}"
            if [ ! -d "${dir_path}" ]; then
                echo "Error: Directory '${dir_path}' does not exist. Please try again."
                dir_path="" # Reset dir_path to force re-prompt
            fi
        done
        echo "${tool_name} ...... OK"
        # Set global variable, e.g., PIGZ_EXE="/home/user/tools/pigz"
        eval "${var_name}='${dir_path}/${tool_name}'"
    fi
}

echo "Checking package dependencies..."
echo ""

# Perform dependency checks for all required tools
check_tool "pigz" "PIGZ_EXE" "/usr/bin"
check_tool "fastp" "FASTP_EXE" "/usr/bin"
check_tool "bbnorm.sh" "BBNORM_SH_EXE" "/path/to/bbtools"
check_tool "redundans.py" "REDUNDANS_PY_EXE" "/path/to/redundans"
check_tool "minimap2" "MINIMAP2_EXE" "/usr/bin"
check_tool "samtools" "SAMTOOLS_EXE" "/usr/bin"
check_tool "runBESST" "RUNBESST_EXE" "/path/to/BESST"
check_tool "GapCloser" "GAPCLOSER_EXE" "/usr/bin"
check_tool "spades.py" "SPADES_PY_EXE" "/usr/bin"
check_tool "seqkit" "SEQKIT_EXE" "/usr/bin"
check_tool "busco" "BUSCO_EXE" "/usr/bin"

echo ""
echo "All package dependencies checked."
echo ""

# --- 4. Helper Function: Split Step Number String ---
# Example: "123" -> "1 2 3"
split_number() {
    local number_string="$1"
    # Use sed to insert a space between each character
    echo "$number_string" | sed 's/\(.\)/\1 /g'
}

# Convert raw steps string to an array
STARTING_STEP_ARRAY=($(split_number "${STARTING_STEPS_RAW}"))

# --- 5. Main Assembly Workflow ---
# Get the absolute path of the species base directory for easier relative path referencing
base_path=$(realpath "${SPECIES_BASE_DIR}")

for num in "${STARTING_STEP_ARRAY[@]}"; do
    echo "--- Running Step ${num} ---"
    current_step_dir="" # Reset current directory variable for each step

    case $num in
        1)
            echo "Step 1: Quality Trimming (fastp)"
            step_name="1-${SPECIES}_fastp"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            # Enter the working directory for the current step
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            "${FASTP_EXE}" -h "${SPECIES}_reads_report.html" \
                -i "${INPUT_R1}" -I "${INPUT_R2}" \
                -o "${prefix1}.fastp.fq.gz" -O "${prefix2}.fastp.fq.gz" \
                --trim_front1=1 --trim_tail1=1 --length_required=20 \
                --dedup --dup_calc_accuracy=6 --correction \
                --trim_poly_g --trim_poly_x \
                --thread "${THREADS}"
            echo ""
            ;;
        2)
            echo "Step 2: Raw Data Normalization (BBtools)"
            step_name="2-${SPECIES}_normalize"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            "${BBNORM_SH_EXE}" in1="../1-${SPECIES}_fastp/${prefix1}.fastp.fq.gz" \
                in2="../1-${SPECIES}_fastp/${prefix2}.fastp.fq.gz" \
                out1="${prefix1}.nor.fq.gz" out2="${prefix2}.nor.fq.gz" \
                target="10" min=2 histcol=2 khist="${SPECIES}_khist.txt" \
                peaks="${SPECIES}_peaks.txt" threads="${THREADS}" -Xmx"${BBTOOLS_MEM_GB}"
            echo ""
            ;;
        3)
            echo "Step 3: Genome Assembly (SPAdes)"
            step_name="3-${SPECIES}_spades"
            # SPAdes creates its output directory within the current directory, so change back to base_path
            cd "${base_path}" || { echo "Error: Could not enter directory ${base_path}"; exit 1; }

            "${SPADES_PY_EXE}" -1 "./2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" \
                -2 "./2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" \
                -o "${step_name}" -k 21,33,55,77 -t "${THREADS}"
            echo ""
            ;;
        4)
            echo "Step 4: Reduction of Heterozygous Contigs (Redundans)"
            step_name="4-${SPECIES}_reduced"
            current_step_dir="${base_path}/${step_name}"
            # Redundans also creates its output directory within the current directory, so change back to base_path
            cd "${base_path}" || { echo "Error: Could not enter directory ${base_path}"; exit 1; }

            "${REDUNDANS_PY_EXE}" -v \
                -f "./3-${SPECIES}_spades/contigs.fasta" \
                -o "${step_name}" \
                -t "${THREADS}" \
                --log "${SPECIES}_${STARTING_STEPS_RAW}_redundans.log" \
                --noscaffolding --nogapclosing --identity 0.7

            # Clean up intermediate files generated by Redundans
            rm -f "${current_step_dir}/contigs"*
            # Move Redundans log file to the log directory
            mv "${SPECIES}_${STARTING_STEPS_RAW}_redundans.log" "${LOG_DIR}/"
            echo ""
            ;;
        5)
            echo "Step 5: Read Mapping (Minimap2 and Samtools)"
            step_name="5-${SPECIES}_map"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            "${MINIMAP2_EXE}" -ax sr \
                "../4-${SPECIES}_reduced/scaffolds.reduced.fa" \
                "../2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" \
                "../2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" \
                -t "${THREADS}" \
                | "${SAMTOOLS_EXE}" sort -@ "${THREADS}" -m "${SAMTOOLS_MEM_PER_THREAD}" -O BAM - -o "${SPECIES}_map.bam"
            "${SAMTOOLS_EXE}" index "${SPECIES}_map.bam"
            echo ""
            ;;
        6)
            echo "Step 6: Scaffolding Assembled Contigs (BESST)"
            step_name="6-${SPECIES}_scaffold"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            "${RUNBESST_EXE}" -c "../4-${SPECIES}_reduced/scaffolds.reduced.fa" \
                -f "../5-${SPECIES}_map/${SPECIES}_map.bam" \
                -o "./" -orientation fr --iter 10000
            echo ""
            ;;
        7)
            echo "Step 7: Gap Filling (GapCloser)"
            step_name="7-${SPECIES}_gapclose"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            # Create gapcloser.config configuration file
            echo "[LIB]" > gapcloser.config
            echo "q1=../2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" >> gapcloser.config
            echo "q2=../2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" >> gapcloser.config

            "${GAPCLOSER_EXE}" -a "../6-${SPECIES}_scaffold/BESST_output/pass1/Scaffolds_pass1.fa" \
                -b gapcloser.config \
                -o "${SPECIES}_scaffolds.gapcloser.fa" \
                -l "${READ_LENGTH}" -t "${THREADS}"

            echo "Generating basic assembly statistics (SeqKit)"
            "${SEQKIT_EXE}" stat -a "${SPECIES}_scaffolds.gapcloser.fa" > assembly.statistics
            echo ""
            ;;
        8)
            echo "Step 8: BUSCO Assessment"
            step_name="8-${SPECIES}_busco"
            current_step_dir="${base_path}/${step_name}"
            mkdir -p "${current_step_dir}"
            cd "${current_step_dir}" || { echo "Error: Could not enter directory ${current_step_dir}"; exit 1; }

            # Check if BUSCO lineage path exists (even if hardcoded, it's good to verify)
            if [ ! -d "${DIR_BUSCO_LINEAGE}" ]; then
                echo "Error: BUSCO lineage directory '${DIR_BUSCO_LINEAGE}' not found. Please ensure the hardcoded path is correct."
                exit 1
            fi

            "${BUSCO_EXE}" -m genome \
                -i "../7-${SPECIES}_gapclose/${SPECIES}_scaffolds.gapcloser.fa" \
                -l "${DIR_BUSCO_LINEAGE}" \
                -c "${THREADS}" \
                -o "${SPECIES}"

            # Copy BUSCO results to the summary directory
            cp -r "./${SPECIES}" "${BUSCO_SUM_DIR}/"
            echo ""
            ;;
        *)
            echo "Warning: Step ${num} is not a valid option. Skipping."
            echo ""
            ;;
    esac
    echo "--- Step ${num} Completed ---"
    echo ""
done

# --- 6. Final Summary ---
end_time=$(date +%s)
cost_time=$(( end_time - start_time ))
echo "--- Script Finished ---"
echo "Total elapsed time: $((cost_time / 3600))h $(((cost_time % 3600) / 60))min $((cost_time % 60))s"
echo "Detailed log file: ${LOG_FILE}"
echo "-----------------------"

exit 0
