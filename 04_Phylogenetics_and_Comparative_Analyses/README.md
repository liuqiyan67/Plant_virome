# 04_phylogenetics_and_comparative_analyses

this directory contains workflows for phylogenetic inference, structural comparison of viral proteins, and bayesian estimation of host transition rates.

## script descriptions

*   **`run_phylogenetic_analysis.sh`**:
    documents the pipeline for multiple sequence alignment (mafft), trimming (trimal), and maximum likelihood tree inference (iq-tree 2). it includes details on the enhanced search parameters used for validation.

*   **`run_dali_structural_comparison.sh`**:
    documents the workflow for all-vs-all protein structure comparison using dalilite, used to infer distant evolutionary relationships for capsid and movement proteins.

*   **`run_empress_cophylogeny.sh`**:
    documents the host-virus cophylogeny analysis using empress. this includes reconciliation analysis and permutation tests for significance.

*   **`run_bayestraits_batch.sh`**:
    a wrapper script that automates the execution of bayestraits v4 for multiple phylogenetic trees. it demonstrates how to run the reversible-jump mcmc (rj-mcmc) analysis for estimating transition rates.

*   **`bayestraits_commands_template.txt`**:
    the template input file containing the specific parameters and priors used for the bayestraits analysis (e.g., iterations, burn-in, rjhp settings).