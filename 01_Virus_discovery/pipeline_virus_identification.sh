#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: pipeline_virus_identification.sh
#
# PURPOSE:
# This script documents the exact bioinformatic commands and parameters used 
# for the discovery and processing of viral contigs from transcriptome data.
#
# NOTE ON REPRODUCIBILITY:
# The commands below represent the actual, functional code executed in our
# HPC environment. File paths and thread counts are specific to our server
# and should be adjusted by the user.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# =============================================================================
# Step 1: De novo transcriptome assembly
# =============================================================================
# Raw paired-end and single-end reads from selected SRA libraries were
# assembled using Trinity v2.14.0 after quality control with Trimmomatic v0.39.

# Example command for a paired-end library:
Trinity --seqType fq --max_memory 150G --left PAIRED_1.fastq --right PAIRED_2.fastq \
        --trimmomatic --output trinity_out_dir --CPU 20 --full_cleanup

# Example command for a single-end library:
Trinity --seqType fq --max_memory 150G --single SINGLE.fastq \
        --trimmomatic --output trinity_out_dir --CPU 20 --full_cleanup


# =============================================================================
# Step 2: Initial contig filtering and processing
# =============================================================================
# Assembled contigs from all sources were combined and subjected to filtering.

# --- A. Removal of rRNA sequences ---
# rRNA sequences were removed by matching the IMG metarRNA database.
# Software: BLAST+ v2.12.0
nohup parallel -j 20 --xapply 'blastn -db /path/to/database/metaSSU_LSU_IMG_IMG/metaSSU_LSU_IMG90_dna.fa -query {1} -out {1}.rrna.out -num_threads 5 -evalue 0.00001 -max_target_seqs 5 -outfmt 6' ::: *.fa &

# --- B. Removal of host contigs ---
# Host-like contigs were identified by BLASTx against the NCBI nr database.
# Software: Diamond v2.0.15
nohup parallel -j 10 --xapply 'diamond blastx --db /path/to/database/nr_db/nr_diamond.dmnd --max-target-seqs 5 --fast --threads 100 --outfmt 6 qseqid sseqid staxids sskingdoms skingdoms sphylums sscinames stitle evalue qstart qend sstart send pident length bitscore mismatch gapopen --query {1} --out {1}.diamond_nr_out' ::: *.fa &

# Note: Contigs with significant hits to non-viral eukaryotic taxa in the nr output 
# were subsequently filtered out using custom parsing scripts.

# --- C. Filter by length ---
# Contigs shorter than 1000 nt were excluded to enrich for more complete genomes.
# Software: seqkit v2.3.0
seqkit seq -m 1000 non_host_contigs.fasta > contigs_gt1000.fasta

# --- D. Cluster at 99% identity ---
# The filtered contigs were clustered to create a set of representative sequences.
# Software: CD-HIT v4.8.1
cd-hit-est -i contigs_gt1000.fasta -o representative_contigs.fasta -c 0.99


# =============================================================================
# Step 3: Iterative virus identification
# =============================================================================
# A five-iteration BLASTp/BLASTx search was performed to identify all potential viral
# RdRp-containing sequences.

# --- Iteration 1: Initial Search ---
# Querying representative contigs against NCBI Virus RdRp reference set.
nohup parallel -j 8 --xapply 'blastx -db ref_rnav_seq.fa -max_target_seqs 100000000 -evalue 0.00001 -num_threads 10 -query {1} -out {1}.orthornavirae.out' ::: *.fas &

# --- Contamination Check (Performed after each iteration) ---
# Hits from each iteration were rigorously checked for contamination before being
# used as queries for the next round.

# Check 1: Ultra-sensitive Diamond BLASTx against nr database to flag non-viral hits.
nohup parallel -j 20 --xapply 'diamond blastx --db /path/to/database/nr_db/nr_diamond.dmnd --max-target-seqs 10 --ultra-sensitive --evalue 0.0000000001 --threads 10 --outfmt 6 qseqid sseqid staxids sskingdoms skingdoms sphylums sscinames stitle evalue qstart qend sstart send pident length bitscore mismatch gapopen --query {1} --out {1}.nr_check.out' ::: *.fan &

# Check 2: InterProScan to identify non-viral protein domains.
# Software: InterProScan v5.65-97.0
nohup parallel -j 6 --xapply '/path/to/interproscan-5.65-97.0/interproscan.sh --cpu 100 -dp --formats TSV,GFF3 --highmem --seqtype n -i {1}' ::: *.nucl &

# Manual Curation Step:
# The outputs from Diamond and InterProScan were manually inspected. 
# - Sequences matching RNA viral similarities or relevant domains were retained.
# - Sequences matching non-viral similarities or domains were removed.
# - Sequences with no similarities were subjected to online BLAST; those matching non-viral genes were removed.


# =============================================================================
# Step 4: Final contig re-assembly
# =============================================================================
# For specific viral groups where fragmented genomes were suspected,
# the identified viral contigs were subjected to a final re-assembly step using CAP3.
# Software: CAP3 (Version date 12/21/07)
cap3 fragment_file.fasta -o 100 -p 98 > cap3.out

echo "--- Pipeline execution complete. ---"
