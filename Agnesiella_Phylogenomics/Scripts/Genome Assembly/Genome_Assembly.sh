#!/bin/bash
# Genome Assembly Pipeline
# 2020.08.04 by ZF
# 2023.07.11 by DSY
# 2024.12.1 by WJJ

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration and Defaults ---
# Default parameters (can be overridden by command-line arguments)
DEFAULT_READ_LENGTH=150
DEFAULT_THREADS=16
DEFAULT_MEMORY_GB=200 # Total memory in GB for tools like SPAdes, BBtools
DEFAULT_BBTOOLS_MEM="80g" # Specific memory setting for BBtools (e.g., 80g)
DEFAULT_SAMTOOLS_SORT_MEM="1706M" # Specific memory setting for Samtools sort (e.g., 1706M)
DEFAULT_FASTP_THREADS=12 # Fastp specific threads, can be linked to THREADS
DEFAULT_BBTOOLS_TARGET_COV=10
DEFAULT_SPADES_K="21,33,55,77"
DEFAULT_REDUNDANS_IDENTITY=0.7
DEFAULT_BESST_ITER=10000
DEFAULT_STEPS="12345678" # Default steps to run

# --- Usage Function ---
usage() {
    echo "Usage: $0 -r1 <read1.fq.gz> -r2 <read2.fq.gz> -s <species_name> -b <busco_lineage_dir> [options]"
    echo "Options:"
    echo "  -r1 <file>             : Input R1 fastq file (required)"
    echo "  -r2 <file>             : Input R2 fastq file (required)"
    echo "  -s <name>              : Species name (required)"
    echo "  -b <dir>               : BUSCO lineage directory (required for step 8)"
    echo "  -l <int>               : Read length (default: $DEFAULT_READ_LENGTH)"
    echo "  -t <int>               : Number of threads (default: $DEFAULT_THREADS)"
    echo "  -m <int>               : Total memory in GB (default: $DEFAULT_MEMORY_GB)"
    echo "  -S <string>            : Steps to run (e.g., '12345678', '345', default: $DEFAULT_STEPS)"
    echo "  -h                     : Show this help message"
    echo ""
    echo "Steps:"
    echo "  1: Quality trimming (fastp)"
    echo "  2: Normalizing coverage (BBtools bbnorm.sh)"
    echo "  3: Genome assembly (SPAdes)"
    echo "  4: Reduction of heterozygous contigs (Redundans)"
    echo "  5: Mapping (minimap2, samtools)"
    echo "  6: Scaffolding (BESST)"
    echo "  7: Gap filling (GapCloser) & Stats (seqkit)"
    echo "  8: BUSCO assessment"
    exit 1
}

# --- Parse Command Line Arguments ---
while getopts "r1:r2:s:b:l:t:m:S:h" opt; do
    case $opt in
        r1) R1_FILE="$OPTARG" ;;
        r2) R2_FILE="$OPTARG" ;;
        s) SPECIES="$OPTARG" ;;
        b) DIR_BUSCO_LINEAGE="$OPTARG" ;;
        l) READ_LENGTH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY_GB="$OPTARG" ;;
        S) STEPS_TO_RUN="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Set defaults if not provided
R1_FILE=${R1_FILE:?"Error: -r1 <read1.fq.gz> is required."}
R2_FILE=${R2_FILE:?"Error: -r2 <read2.fq.gz> is required."}
SPECIES=${SPECIES:?"Error: -s <species_name> is required."}
DIR_BUSCO_LINEAGE=${DIR_BUSCO_LINEAGE:?"Error: -b <busco_lineage_dir> is required for step 8."} # Make required if step 8 is included
READ_LENGTH=${READ_LENGTH:-$DEFAULT_READ_LENGTH}
THREADS=${THREADS:-$DEFAULT_THREADS}
MEMORY_GB=${MEMORY_GB:-$DEFAULT_MEMORY_GB}
STEPS_TO_RUN=${STEPS_TO_RUN:-$DEFAULT_STEPS}

# Specific memory settings derived from total memory or defaults
BBTOOLS_MEM="${MEMORY_GB}g" # Use total memory for BBtools Xmx
SAMTOOLS_SORT_MEM="$DEFAULT_SAMTOOLS_SORT_MEM" # Keep default or make configurable
FASTP_THREADS="$DEFAULT_FASTP_THREADS" # Keep default or link to THREADS

# --- Dependency Check Function ---
# Checks if a tool is in PATH or prompts user for directory
# Args: $1=tool_name, $2=description, $3=variable_name_for_path
check_tool() {
    local tool_name="$1"
    local description="$2"
    local path_var_name="$3"
    local tool_path=""

    echo "Checking for $description ($tool_name)......"

    # Check if tool is in PATH
    if tool_path=$(which "$tool_name"); then
        echo "$description ...... OK ($tool_path)"
        eval "export $path_var_name='$tool_path'"
        return 0
    else
        echo "$description not found in PATH."
        # Prompt user for directory until valid executable is found
        while true; do
            read -p "Please input the installation directory for $description (absolute path, e.g. /usr/bin): " user_dir
            if [ -x "$user_dir/$tool_name" ]; then
                tool_path="$user_dir/$tool_name"
                echo "$description ...... OK ($tool_path)"
                eval "export $path_var_name='$tool_path'"
                return 0
            else
                echo "Error: $user_dir/$tool_name is not an executable file. Please try again."
            fi
        done
    fi
}

# --- Check Dependencies ---
echo "Checking package dependencies......"
check_tool pigz "pigz" PATH_PIGZ
check_tool fastp "Fastp" PATH_FASTP
check_tool bbnorm.sh "BBtools bbnorm.sh" PATH_BBNORM
check_tool redundans.py "Redundans" PATH_REDUNDANS
check_tool minimap2 "Minimap2" PATH_MINIMAP2
check_tool samtools "Samtools" PATH_SAMTOOLS
check_tool runBESST "BESST runBESST" PATH_BESST
check_tool GapCloser "GapCloser" PATH_GAPCLOSER
check_tool seqkit "SeqKit" PATH_SEQKIT # Added seqkit check
check_tool busco "BUSCO" PATH_BUSCO # Added busco check

echo "" # Add a newline after checks

# --- Setup Directories and Logging ---
start_time=$(date +%s)

# Create main directories
mkdir -p 0-log 1-busco_sum "${SPECIES}"

# Define base path for species-specific output
BASE_PATH=$(realpath "${SPECIES}")

# Define log file path
LOG_FILE="./0-log/${SPECIES}_${STEPS_TO_RUN}_assemble.log"

# Redirect all subsequent output (stdout and stderr) to the log file
# This should be done *after* initial setup and dependency checks if you want
# those initial messages to go to stdout as well. Or do it earlier if you want
# everything in the log. Let's put everything in the log from here.
exec &>> "$LOG_FILE"

echo "--- Pipeline Start ---"
echo "Start Time: $(date)"
echo "Input R1: $(realpath "$R1_FILE")"
echo "Input R2: $(realpath "$R2_FILE")"
echo "Species: $SPECIES"
echo "BUSCO Lineage Dir: $(realpath "$DIR_BUSCO_LINEAGE")"
echo "Read Length: $READ_LENGTH"
echo "Threads: $THREADS"
echo "Total Memory (GB): $MEMORY_GB"
echo "Steps to Run: $STEPS_TO_RUN"
echo "Base Output Path: $BASE_PATH"
echo ""

# --- Prepare File Prefixes ---
filename1=$(basename "$R1_FILE")
filename2=$(basename "$R2_FILE")
prefix1="${filename1%%.*}"
prefix2="${filename2%%.*}"

# --- Execute Steps ---
# Split steps string into an array
IFS='' read -r -a STARTING_STEP <<< "$(echo "$STEPS_TO_RUN" | sed 's/./& /g')"

for num in "${STARTING_STEP[@]}"
do
    echo "--- Running Step $num ---"
    step_start_time=$(date +%s)

    case $num in
		1)
			echo "Step 1: Quality trimming (fastp)"
			STEP_DIR="${BASE_PATH}/1-${SPECIES}_fastp"
			mkdir -p "$STEP_DIR"

			"$PATH_FASTP" \
                -h "${STEP_DIR}/${SPECIES}_reads_report.html" \
                -i "$R1_FILE" \
                -I "$R2_FILE" \
                -o "${STEP_DIR}/${prefix1}.fastp.fq.gz" \
                -O "${STEP_DIR}/${prefix2}.fastp.fq.gz" \
                --trim_front1=1 \
                --trim_tail1=1 \
                --length_required=20 \
                --dedup \
                --dup_calc_accuracy=6 \
                --correction \
                --trim_poly_g \
                --trim_poly_x \
                --thread "$FASTP_THREADS"
            ;;

		2)
			echo "Step 2: Normalizing coverage of raw data (BBtools bbnorm.sh)"
			STEP_DIR="${BASE_PATH}/2-${SPECIES}_normalize"
			mkdir -p "$STEP_DIR"

            "$PATH_BBNORM" \
                in1="${BASE_PATH}/1-${SPECIES}_fastp/${prefix1}.fastp.fq.gz" \
                in2="${BASE_PATH}/1-${SPECIES}_fastp/${prefix2}.fastp.fq.gz" \
                out1="${STEP_DIR}/${prefix1}.nor.fq.gz" \
                out2="${STEP_DIR}/${prefix2}.nor.fq.gz" \
                target="$DEFAULT_BBTOOLS_TARGET_COV" \
                min=2 \
                histcol=2 \
                khist="${STEP_DIR}/${SPECIES}_khist.txt" \
                peaks="${STEP_DIR}/${SPECIES}_peaks.txt" \
                threads="$THREADS" \
                -Xmx"$BBTOOLS_MEM"
            ;;

		3)
			echo "Step 3: Genome assembly (SPAdes)"
			STEP_DIR="${BASE_PATH}/3-${SPECIES}_spades"
			# SPAdes creates its own directory, just ensure BASE_PATH exists
			mkdir -p "$BASE_PATH"

			"$PATH_SPADES" \
                -1 "${BASE_PATH}/2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" \
                -2 "${BASE_PATH}/2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" \
                -o "$STEP_DIR" \
                -k "$DEFAULT_SPADES_K" \
                -t "$THREADS" \
                --memory "$MEMORY_GB" # Use total memory variable
            ;;

		4)
			echo "Step 4: Reduction of heterozygous contigs (Redundans)"
			STEP_DIR="${BASE_PATH}/4-${SPECIES}_reduced"
			mkdir -p "$STEP_DIR"

			"$PATH_REDUNDANS" \
                -v \
                -f "${BASE_PATH}/3-${SPECIES}_spades/contigs.fasta" \
                -o "$STEP_DIR" \
                -t "$THREADS" \
                --log "${STEP_DIR}/${SPECIES}_${STEPS_TO_RUN}_redundans.log" \
                --noscaffolding \
                --nogapclosing \
                --identity "$DEFAULT_REDUNDANS_IDENTITY"

            # Clean up intermediate files if desired
			# rm "${STEP_DIR}/contigs*" # Uncomment if you want to remove them

            # Move redundans log to main log directory
            mv "${STEP_DIR}/${SPECIES}_${STEPS_TO_RUN}_redundans.log" "./0-log/"
			;;

		5)
			echo "Step 5: Mapping (minimap2, samtools)"
			STEP_DIR="${BASE_PATH}/5-${SPECIES}_map"
			mkdir -p "$STEP_DIR"

			"$PATH_MINIMAP2" \
                -ax sr \
                "${BASE_PATH}/4-${SPECIES}_reduced/scaffolds.reduced.fa" \
                "${BASE_PATH}/2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" \
                "${BASE_PATH}/2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" \
                -t "$THREADS" \
            | "$PATH_SAMTOOLS" sort \
                -@ "$THREADS" \
                -m "$SAMTOOLS_SORT_MEM" \
                -O BAM \
                - \
                -o "${STEP_DIR}/${SPECIES}_map.bam"

            "$PATH_SAMTOOLS" index "${STEP_DIR}/${SPECIES}_map.bam"
            ;;

		6)
			echo "Step 6: Scaffolding assembled contigs (BESST)"
			STEP_DIR="${BASE_PATH}/6-${SPECIES}_scaffold"
			mkdir -p "$STEP_DIR"

			"$PATH_BESST" \
                -c "${BASE_PATH}/4-${SPECIES}_reduced/scaffolds.reduced.fa" \
                -f "${BASE_PATH}/5-${SPECIES}_map/${SPECIES}_map.bam" \
                -o "$STEP_DIR" \
                -orientation fr \
                --iter "$DEFAULT_BESST_ITER"
            ;;

		7)
			echo "Step 7: Gap filling (GapCloser) & Generate basic assembly statistics (seqkit)"
			STEP_DIR="${BASE_PATH}/7-${SPECIES}_gapclose"
			mkdir -p "$STEP_DIR"

            # Create gapcloser config file
            GAPCLOSER_CONFIG="${STEP_DIR}/gapcloser.config"
			echo "[LIB]" > "$GAPCLOSER_CONFIG"
			echo "q1=${BASE_PATH}/2-${SPECIES}_normalize/${prefix1}.nor.fq.gz" >> "$GAPCLOSER_CONFIG"
			echo "q2=${BASE_PATH}/2-${SPECIES}_normalize/${prefix2}.nor.fq.gz" >> "$GAPCLOSER_CONFIG"

            # Run GapCloser
			"$PATH_GAPCLOSER" \
                -a "${BASE_PATH}/6-${SPECIES}_scaffold/BESST_output/pass1/Scaffolds_pass1.fa" \
                -b "$GAPCLOSER_CONFIG" \
                -o "${STEP_DIR}/${SPECIES}_scaffolds.gapcloser.fa" \
                -l "$READ_LENGTH" \
                -t "$THREADS"

			echo "Generating basic assembly statistics using seqkit..."
			"$PATH_SEQKIT" stat -a "${STEP_DIR}/${SPECIES}_scaffolds.gapcloser.fa" > "${STEP_DIR}/assembly.statistics"
            ;;

		8)
			echo "Step 8: BUSCO assessment"
			STEP_DIR="${BASE_PATH}/8-${SPECIES}_busco"
			mkdir -p "$STEP_DIR"

			"$PATH_BUSCO" \
                -m genome \
                -i "${BASE_PATH}/7-${SPECIES}_gapclose/${SPECIES}_scaffolds.gapcloser.fa" \
                -l "$DIR_BUSCO_LINEAGE" \
                -c "$THREADS" \
                -o "$SPECIES" \
                --out_dir "$STEP_DIR" # Use --out_dir to specify output location

            # Copy BUSCO results summary to the main summary directory
			cp -r "${STEP_DIR}/run_${SPECIES}" "../1-busco_sum/"
			;;

		*)
            echo "Warning: Unknown step number '$num'. Skipping."
            ;;
    esac

    step_end_time=$(date +%s)
    step_cost_time=$[ $step_end_time-$step_start_time ]
    echo "--- Step $num finished in $(($step_cost_time / 3600))h $((step_cost_time % 3600 / 60))min $(($step_cost_time % 60))s ---"
    echo "" # Add newline after each step summary
done

# --- Final Timing ---
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "--- Pipeline Finished ---"
echo "End Time: $(date)"
echo "Total elapsed time: $(($cost_time / 3600))h $((cost_time % 3600 / 60))min $(($cost_time % 60))s"
echo "Log file: $(realpath "$LOG_FILE")"
echo "--- End of Log ---"
