#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: calculate_virus_rpm.sh
#
# PURPOSE:
# This script documents the workflow for calculating viral abundance (Reads Per
# Million, RPM) by mapping cleaned RNA-seq reads to viral genomes.
#
# The workflow consists of:
# 1. Identifying the target virus for each SRA library.
# 2. Mapping reads to the viral genome using Bowtie2.
# 3. Converting SAM output to sorted BAM files using Samtools.
# 4. Calculating mapping statistics (mapped reads count, total bases).
#
# INPUT:
#   - Cleaned FASTQ files (e.g., from Trimmomatic)
#   - Viral genome indices (built with bowtie2-build)
#   - A mapping file linking SRA IDs to Virus IDs
#
# OUTPUT:
#   - Sorted BAM files
#   - A summary table of viral abundance statistics
#
# USAGE:
# This script is for documentation. To run, replace placeholder paths.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Step 0: Define Variables ---
# Directory containing trimmed FASTQ files
TRIM_DIR="/path/to/trimmed_fastq"
# Directory containing Bowtie2 indices for viral genomes
INDEX_DIR="/path/to/viral_indices"
# Output directory for BAM files
BAM_DIR="/path/to/output_bam"
# Tab-delimited file linking SRA ID (col 1) to Virus ID (col 2)
MAP_FILE="/path/to/sra_virus_map.txt"

mkdir -p "${BAM_DIR}"


# =============================================================================
# Step 1: Mapping Reads to Viral Genomes
# =============================================================================
echo "Step 1: Mapping reads using Bowtie2..."

# Loop through the mapping file to process each SRA-Virus pair
while read -r SRA_ID VIRUS_ID; do
    
    FASTQ_FILE="${TRIM_DIR}/${SRA_ID}_trimmed.fq"
    INDEX_PREFIX="${INDEX_DIR}/${VIRUS_ID}"
    SAM_OUTPUT="${BAM_DIR}/${SRA_ID}_${VIRUS_ID}.sam"
    BAM_OUTPUT="${BAM_DIR}/${SRA_ID}_${VIRUS_ID}.bam"
    
    # Check if input files exist
    if [[ -f "${FASTQ_FILE}" && -f "${INDEX_PREFIX}.1.bt2" ]]; then
        
        echo "Processing: ${SRA_ID} -> ${VIRUS_ID}"
        
        # 1. Run Bowtie2 mapping
        # --local: Perform local alignment (useful for potentially divergent viral reads)
        # -L 20 -N 1: Sensitivity settings
        bowtie2 -x "${INDEX_PREFIX}" \
                -U "${FASTQ_FILE}" \
                --local -L 20 -N 1 \
                --threads 8 \
                -S "${SAM_OUTPUT}" 2> "${BAM_DIR}/${SRA_ID}.log"
        
        # 2. Convert SAM to sorted BAM using Samtools
        # -F 4: Exclude unmapped reads to save space
        samtools view -bS -F 4 "${SAM_OUTPUT}" | \
        samtools sort -o "${BAM_OUTPUT}"
        
        # 3. Index the BAM file
        samtools index "${BAM_OUTPUT}"
        
        # Remove the large intermediate SAM file
        rm "${SAM_OUTPUT}"
        
    else
        echo "Warning: Input file or index missing for ${SRA_ID}"
    fi

done < "${MAP_FILE}"


# =============================================================================
# Step 2: Calculating Abundance Statistics
# =============================================================================
echo "Step 2: Calculating mapped reads and abundance statistics..."

OUTPUT_STATS="${BAM_DIR}/virus_abundance.txt"
echo -e "SRA_ID\tVirus_ID\tMapped_Reads\tMean_Read_Length\tTotal_Mapped_Bases" > "${OUTPUT_STATS}"

for BAM_FILE in ${BAM_DIR}/*.bam; do
    
    # Extract IDs from filename
    BASENAME=$(basename "${BAM_FILE}" .bam)
    # Assumes filename format: SRAID_VirusID.bam
    SRA_ID=$(echo "${BASENAME}" | cut -d_ -f1)
    VIRUS_ID=$(echo "${BASENAME}" | cut -d_ -f2-)
    
    # Count mapped reads
    MAPPED_READS=$(samtools view -c "${BAM_FILE}")
    
    # Calculate total mapped bases and mean read length
    # Using awk to sum the length of the sequence field (column 10)
    STATS=$(samtools view "${BAM_FILE}" | awk '{len+=length($10); cnt++} END {if (cnt>0) print len/cnt, len; else print 0, 0}')
    MEAN_LEN=$(echo "${STATS}" | cut -d' ' -f1)
    TOTAL_BASES=$(echo "${STATS}" | cut -d' ' -f2)
    
    echo -e "${SRA_ID}\t${VIRUS_ID}\t${MAPPED_READS}\t${MEAN_LEN}\t${TOTAL_BASES}" >> "${OUTPUT_STATS}"

done

echo "--- Abundance calculation workflow complete! ---"