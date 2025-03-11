#!/bin/bash
# Enhanced GATK Pipeline Script v2024.12.1
# Author: WJJ

set -euo pipefail

# Global timer initialization
MAIN_START=$(date +%s)

# Parameter assignment
FASTQ_R1="$1"
FASTQ_R2="$2"
SPECIES_NAME="$3"
STEPS="$4"

# Initialize logging
LOG_DIR="0-log"
mkdir -p "${LOG_DIR}"
LOG="${LOG_DIR}/${SPECIES_NAME}_${STEPS}_pipeline.log"
echo "Pipeline started: $(date)" > "${LOG}"
echo "species_name=${SPECIES_NAME}" >> "${LOG}"

# Helper functions
split_steps() {
    grep -o . <<< "$1" | tr '\n' ' '
}

timer() {
    local start=$1
    local end=$(date +%s)
    echo "$((end - start)) seconds"
}

# Tool check
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "ERROR: $1 not found in PATH" | tee -a "${LOG}"
        exit 1
    fi
}

# Check essential tools
check_tool bwa-mem2
check_tool samtools
check_tool gatk
check_tool bcftools
check_tool vcftools

# Reference genome path
REF="/path/to/reference/genome.fasta"

# Create workspace
mkdir -p "1-bams" "2-g_vcfs" "3-vcfs" "4-snps" "5-snps_filted" "6-snps_further_filted"

# Step executor
execute_step() {
    local step_num=$1
    local step_start=$(date +%s)
    
    case $step_num in
    1)
        # Alignment and Sorting
        echo "[$(date)] STEP 1 START: Alignment" | tee -a "${LOG}"
        cd 1-bams

        bwa-mem2 mem "${REF}" "${FASTQ_R1}" "${FASTQ_R2}" \
            -R "@RG\tID:${SPECIES_NAME}\tSM:${SPECIES_NAME}" \
            -M -t 24 | samtools sort -@ 12 -m 1706M -O BAM \
            -o "${SPECIES_NAME}_sorted.bam"

        cd ..
        ;;
    2)
        # Duplicate Marking
        echo "[$(date)] STEP 2 START: Mark Duplicates" | tee -a "${LOG}"
        cd 1-bams

        gatk --java-options "-Xmx80G -XX:ParallelGCThreads=16" MarkDuplicates \
            -I "${SPECIES_NAME}_sorted.bam" \
            -O "${SPECIES_NAME}_markdup.bam" \
            -M "${SPECIES_NAME}_markdup_metrics.txt" \
            --REMOVE_DUPLICATES true

        gatk BuildBamIndex -I "${SPECIES_NAME}_markdup.bam"
        cd ..
        ;;
    3)
        # gVCF Generation
        echo "[$(date)] STEP 3 START: gVCF Creation" | tee -a "${LOG}"
        cd 2-g_vcfs

        gatk HaplotypeCaller \
            -R "${REF}" \
            -I ../1-bams/"${SPECIES_NAME}_markdup.bam" \
            --max-alternate-alleles 3 \
            --sample-ploidy 2 \
            -ERC GVCF \
            -O "${SPECIES_NAME}.g.vcf.gz"

        tabix -p vcf "${SPECIES_NAME}.g.vcf.gz"
        cd ..
        ;;
    4)
        # Joint Genotyping
        echo "[$(date)] STEP 4 START: Genotyping" | tee -a "${LOG}"
        cd 3-vcfs

        gatk --java-options "-Xmx4g" GenotypeGVCFs \
            -R "${REF}" \
            -V ../2-g_vcfs/"${SPECIES_NAME}.g.vcf.gz" \
            -O "${SPECIES_NAME}.vcf.gz"

        tabix -p vcf "${SPECIES_NAME}.vcf.gz"
        cd ..
        ;;
    5)
        # SNP Selection
        echo "[$(date)] STEP 5 START: SNP Selection" | tee -a "${LOG}"
        cd 4-snps

        gatk --java-options "-Xmx4G" SelectVariants \
            -R "${REF}" \
            -V ../3-vcfs/"${SPECIES_NAME}.vcf.gz" \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -O "${SPECIES_NAME}_BIALLELIC_SNP.vcf.gz"

        cd ..
        ;;
    6)
        # SNP Filtering
        echo "[$(date)] STEP 6 START: SNP Filtering" | tee -a "${LOG}"
        cd 5-snps_filted

        gatk VariantFiltration \
            -V "../4-snps/${SPECIES_NAME}_BIALLELIC_SNP.vcf.gz" \
            --filter-expression "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --filter-name "my_filters" \
            -O "${SPECIES_NAME}_BIALLELIC_SNP_anno.vcf.gz"

        gatk --java-options "-Xmx4G" SelectVariants \
            -R "${REF}" \
            -V "${SPECIES_NAME}_BIALLELIC_SNP_anno.vcf.gz" \
            --exclude-filtered \
            -O "${SPECIES_NAME}_BIALLELIC_SNP_PASS.vcf.gz"

        tabix -p vcf "${SPECIES_NAME}_BIALLELIC_SNP_PASS.vcf.gz"
        cd ..
        ;;
    7)
        # VCF Merging
        echo "[$(date)] STEP 7 START: VCF Merging" | tee -a "${LOG}"
        cd 5-snps_filted

        bcftools merge *_BIALLELIC_SNP_PASS.vcf.gz \
            -o "All_BIALLELIC_SNP_PASS.vcf.gz"

        tabix -p vcf "All_BIALLELIC_SNP_PASS.vcf.gz"
        cd ..
        ;;
    8)
        # Advanced Filtering
        echo "[$(date)] STEP 8 START: Advanced Filtering" | tee -a "${LOG}"
        cd 6-snps_further_filted

        vcftools --gzvcf "../5-snps_filted/All_BIALLELIC_SNP_PASS.vcf.gz" \
            --recode --recode-INFO-all --stdout \
            --max-missing 0.8 --minDP 4 --maxDP 100 \
            --minQ 30 --min-alleles 2 --max-alleles 2 | bgzip -c > "All_BIALLELIC_SNP_PASS_0.8_30.vcf.gz"

        bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.5' \
            -m2 -M2 -o "All_BIALLELIC_SNP_PASS_0.8_30_0.5.vcf.gz" \
            -O z ./"All_BIALLELIC_SNP_PASS_0.8_30.vcf.gz"
        cd ..
        ;;
    *)
        echo "WARNING: Ignored invalid step ${step_num}" | tee -a "${LOG}"
        return
        ;;
    esac

    local duration=$(timer "$step_start")
    echo "[$(date)] STEP ${step_num} COMPLETED: Time ${duration}" | tee -a "${LOG}"
}

# Execute requested steps
for step in $(split_steps "${STEPS}"); do
    execute_step "$step"
done

# Final timing statistics
TOTAL_DURATION=$(timer "$MAIN_START")
echo "PROCESSING SUMMARY:" | tee -a "${LOG}"
echo "Total pipeline time: ${TOTAL_DURATION}" | tee -a "${LOG}"
echo "Execution completed at: $(date)" | tee -a "${LOG}"
