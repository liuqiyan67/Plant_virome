#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: siRNA_analysis_mapping.sh
#
# PURPOSE:
# This script documents the workflow for mapping small RNA (sRNA) reads to
# viral genomes and calculating coverage statistics to validate host association.
#
# The workflow consists of:
# 1. Mapping cleaned sRNA reads to a viral genome index using Bowtie2.
# 2. Sorting and indexing the resulting BAM files.
# 3. Calculating detailed coverage and depth statistics for each viral contig.
#
# INPUT:
#   - Cleaned sRNA FASTQ files (from the preprocessing step)
#   - Bowtie2 index of the viral genomes
#
# OUTPUT:
#   - Sorted BAM files
#   - A summary table (virus_reference_stats.txt) containing mapped read counts,
#     genome coverage (%), and mean depth.
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Directory containing cleaned sRNA FASTQ files
FASTQ_DIR="/path/to/clean_fastq"
# Bowtie2 index basename for the viral genomes
INDEX_BASE="/path/to/viral_index/viral_db"
# Output directory for alignment results
ALIGN_DIR="/path/to/alignment_output"
# Output file for the final statistics table
STATS_FILE="${ALIGN_DIR}/virus_reference_stats.txt"

mkdir -p "${ALIGN_DIR}"


# =============================================================================
# Step 1: Bowtie2 Alignment
# =============================================================================
echo "Step 1: Mapping sRNA reads to viral genomes..."

# Loop through all cleaned FASTQ files
for FASTQ in "${FASTQ_DIR}"/*_clean.fq; do
    
    BASENAME=$(basename "${FASTQ}" _clean.fq)
    BAM_OUTPUT="${ALIGN_DIR}/${BASENAME}.bam"
    LOG_FILE="${ALIGN_DIR}/${BASENAME}.log"
    
    echo "Aligning: ${BASENAME}"
    
    # Run Bowtie2
    # --very-sensitive-local: Accurate alignment for short reads
    # -N 1 -L 15: Allow 1 mismatch in seed, shorter seed length (crucial for sRNA)
    # -k 50: Report up to 50 alignments per read (useful for multi-mapping analysis)
    bowtie2 --very-sensitive-local -x "${INDEX_BASE}" -U "${FASTQ}" \
            -N 1 -L 15 -k 50 \
            --threads 8 \
            2> "${LOG_FILE}" | \
    samtools view -bS - | \
    samtools sort -o "${BAM_OUTPUT}"
    
    # Index the sorted BAM file
    samtools index "${BAM_OUTPUT}"
    
done


# =============================================================================
# Step 2: Calculate Coverage and Abundance Statistics
# =============================================================================
echo "Step 2: Calculating coverage statistics..."

# Initialize the output table with headers
echo -e "Sample_ID\tReference_ID\tSeq_Length\tMapped_Reads\tCoverage_Rate\tMean_Depth" > "${STATS_FILE}"

# Loop through all generated BAM files
for BAM in "${ALIGN_DIR}"/*.bam; do
    
    SAMPLE_ID=$(basename "${BAM}" .bam)
    echo "Processing stats for: ${SAMPLE_ID}"
    
    # Get per-reference statistics using samtools idxstats
    # idxstats output format: ref_name seq_length mapped_reads unmapped_reads
    samtools idxstats "${BAM}" | while read -r REF_ID REF_LEN MAPPED UNMAPPED; do
        
        # Skip the unmapped entry ('*') or empty lines
        if [[ "${REF_ID}" == "*" || -z "${REF_ID}" ]]; then continue; fi
        
        # Only process references that have mapped reads
        if [[ "${MAPPED}" -gt 0 ]]; then
            
            # Calculate depth and coverage using samtools depth
            # This calculates the depth at each position
            samtools depth -r "${REF_ID}" "${BAM}" > "${ALIGN_DIR}/temp.depth"
            
            # Calculate Covered Bases (count of positions with depth > 0)
            COVERED_BASES=$(awk '$3 > 0 {count++} END {print count+0}' "${ALIGN_DIR}/temp.depth")
            
            # Calculate Total Depth (sum of depth at all positions)
            TOTAL_DEPTH=$(awk '{sum+=$3} END {print sum+0}' "${ALIGN_DIR}/temp.depth")
            
            # Compute Coverage Rate (Covered Bases / Total Length)
            COVERAGE_RATE=$(awk -v cov="${COVERED_BASES}" -v len="${REF_LEN}" 'BEGIN {printf "%.6f", cov/len}')
            
            # Compute Mean Depth (Total Depth / Total Length)
            MEAN_DEPTH=$(awk -v sum="${TOTAL_DEPTH}" -v len="${REF_LEN}" 'BEGIN {printf "%.2f", sum/len}')
            
            # Append results to the statistics file
            echo -e "${SAMPLE_ID}\t${REF_ID}\t${REF_LEN}\t${MAPPED}\t${COVERAGE_RATE}\t${MEAN_DEPTH}" >> "${STATS_FILE}"
            
            # Clean up temp file
            rm "${ALIGN_DIR}/temp.depth"
        fi
    done
done

echo "--- Mapping and analysis complete. Results saved to ${STATS_FILE} ---"