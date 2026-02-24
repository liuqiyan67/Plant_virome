# Analysis pipeline for "The Origins and Evolution of the RNA Virome in Land Flora"

This repository contains the custom analysis scripts and bioinformatic workflows used in the manuscript.

## Purpose of this repository

The scripts provided here are intended to document the computational methods used in our study. They serve as a transparent record of the exact parameters, logic, and sequence of operations performed.

**Please note:** these scripts are not designed as a standalone software package. They are not intended to be executed "as-is" in any arbitrary environment. To run them, users would need to:
1.  Install the required software dependencies (listed below).
2.  Adjust file paths and system-specific configurations (e.g., path to databases, number of threads) to match their own computing environment.

## Workflow overview

The analysis is organized into four main stages, corresponding to the numbered directories:

*   **01_virus_discovery**: pipeline for de novo assembly and identification of viral contigs from transcriptome data.
*   **02_quality_control**: the multi-tiered framework for assessing host veracity, including library purity checks and siRNA validation.
*   **03_diversity_and_network_analyses**: scripts for clustering, rarefaction analysis, and constructing sequence similarity networks.
*   **04_phylogenetics_and_comparative_analyses**: workflows for phylogenetic tree inference, structural comparisons, and bayesian transition rate estimation.

## Software Dependencies

The following software tools and specific versions were used in this study. While these are also listed in the manuscript's Key Resources Table, they are provided here for immediate reference to facilitate environment setup and reproducibility.

**Quality Control & De Novo Assembly**
* Trimmomatic v0.39
* Trinityrnaseq v2.14.0
* CAP3

**Sequence Processing & Alignment**
* SeqKit v2.0.0
* Bowtie2
* BWA
* MAFFT
* trimAl v1.2

**Virus Identification & Similarity Search**
* Diamond v2.0.15
* NCBI BLAST+ suite v2.12.0+
* HMMER v3.3

**Clustering & Taxonomy**
* CD-HIT v4.8.1
* MMseqs2 v13.45111
* MCL v14-137

**Phylogenetics & Evolutionary Analysis**
* IQ-TREE v2
* FastTree
* BayesTraits V4

**Protein Structure Prediction & Comparison**
* AlphaFold3
* Phyre2
* DALI

**Data Processing, Statistics & Visualization**
* GNU Parallel v1
* Python 3 (pandas)
* R (tidyverse, ggplot2, igraph, iNEXT, vegan)
* iTOL
* Cirit

## Data availability

All raw data, intermediate files (such as alignments and tree files), and final result tables associated with this study are publicly available at our figshare repository: **https://doi.org/10.6084/m9.figshare.30688322.v1**

## Contact

For questions regarding the code or analysis, please contact Guan-Zhu Han at guanzhu@njnu.edu.cn






