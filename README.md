# Analysis pipeline for "the origins and evolution of the RNA virome in land flora"

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

## Software dependencies

The following major software tools were used in this study. Please refer to the methods section of the manuscript for specific version numbers.

*   **data processing**: python 3, r, pandas, tidyverse
*   **assembly & identification**: trinity, diamond, mmseqs2
*   **phylogenetics**: iq-tree 2, mafft, trimal
*   **comparative analysis**: bayestraits v4, dalilite
*   **visualization**: ggplot2, igraph, itol

## Data availability

All raw data, intermediate files (such as alignments and tree files), and final result tables associated with this study are publicly available at our figshare repository: https://figshare.com/s/404e77121795e94cce7a

## Contact

For questions regarding the code or analysis, please contact Guan-Zhu Han at guanzhu@njnu.edu.cn

