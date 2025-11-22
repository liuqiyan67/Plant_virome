#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: pipeline_virus_identification.sh
#
# PURPOSE:
# This script documents the key bioinformatic steps for the discovery
# and processing of viral contigs from transcriptome data, as described in
# our manuscript. It is intended as a guide to the methodology, not as a

# fully automated pipeline.
#
# Each step represents a conceptual stage of the analysis, with example
# commands provided.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo "--- This is a documentation script, not intended for direct execution without modification. ---"
echo "--- Please see the README and individual comments for context. ---"

# =============================================================================
# Step 1: De novo transcriptome assembly
# =============================================================================
# Raw paired-end and single-end reads from selected SRA libraries were
# assembled using Trinity after quality control with Trimmomatic.

# Example command for a paired-end library:
Trinity --seqType fq --max_memory 150G --left PAIRED_1.fastq --right PAIRED_2.fastq \
        --trimmomatic --output trinity_out_dir --CPU 20 --full_cleanup

# Example command for a single-end library:
Trinity --seqType fq --max_memory 150G --single SINGLE.fastq \
        --trimmomatic --output trinity_out_dir --CPU 20 --full_cleanup


# =============================================================================
# Step 2: Initial contig filtering and processing
# =============================================================================
# Assembled contigs from all sources were combined and subjected to a
# series of filtering steps to remove non-viral and redundant sequences.

# Assumed starting file: all_assembled_contigs.fasta

# --- A. Removal of host and rRNA contigs (Conceptual) ---
# Host-like contigs were identified by BLASTx against the NCBI nr database and removed.
# (Example command not shown for brevity, as it's a standard procedure).

# --- B. Filter by length ---
# Contigs shorter than 1000 nt were excluded to enrich for more complete genomes.
echo "Step 2B: Filtering contigs by length..."
seqkit seq -m 1000 non_host_contigs.fasta > contigs_gt1000.fasta

# --- C. Cluster at 99% identity ---
# The filtered contigs were clustered to create a set of representative sequences.
echo "Step 2C: Clustering contigs at 99% identity..."
cd-hit-est -i contigs_gt1000.fasta -o representative_contigs.fasta -c 0.99


# =============================================================================
# Step 3: Iterative virus identification
# =============================================================================
# A five-iteration BLASTp search was performed to identify all potential viral
# RdRp-containing sequences from the translated representative contigs.

# First, open reading frames (ORFs) were predicted from the representative contigs.
# Example using ORFfinder:
orffinder -in representative_contigs.fasta -ml 450 -out representative_proteins.fasta

# The iterative process is described conceptually below.
echo "Step 3: Performing iterative BLASTp (conceptual)..."

# --- Contamination check within each iteration ---
# Hits from each iteration were carefully checked for contamination before being
# used as queries for the next round. This involved two checks:
echo "  -> Contamination check for iteration hits..."

# Check 1: BLASTp against nr database to flag non-viral hits.
# Check 2: InterProScan to identify non-viral protein domains.

# (These outputs were manually inspected to filter out any sequences showing
# strong hits to non-viral eukaryotic protein families.)


# =============================================================================
# Step 4: Final contig re-assembly
# =============================================================================
# For specific viral groups where fragmented genomes were suspected,
# the identified viral contigs were subjected to a final re-assembly step using CAP3.
echo "Step 4: Performing CAP3 re-assembly on specific subsets..."
cap3 your_fragment_file.fasta -o 100 -i 98

echo "--- Documentation of the main discovery pipeline is complete. ---"