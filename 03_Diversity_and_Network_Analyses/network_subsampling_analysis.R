# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: network_subsampling_analysis.R
#
# PURPOSE:
# This script performs a stratified subsampling experiment to test the robustness
# of network modularity against uneven sampling of different viral phyla.
#
# It creates balanced sub-networks by sampling an equal number of nodes from
# each major phylum, recalculates modularity, and compares the distribution
# of these values to the modularity of the original, full network.
#
# INPUT:
#   - Edge list: A tab-delimited file (Source, Target)
#   - Node metadata: An Excel/CSV file with node IDs and Phylum information
#
# OUTPUT:
#   - A CSV file containing modularity scores from all replicates.
#   - A PDF plot visualizing the comparison.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. LOAD LIBRARIES ----
# ----------------------
if (!require("igraph")) install.packages("igraph")
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggdist")) install.packages("ggdist")

library(igraph)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)

# 2. SETUP: PARAMETERS AND FILE PATHS ----
# ----------------------------------------
# !!! PLEASE UPDATE THESE PATHS !!!
# For reproducibility, set a working directory or use relative paths
# setwd("path/to/your/data")

edge_file <- "all_plantsrelated_rdrpcores_mmseqclu9.fa.mmseqclu5.fa.e-10.blastpout.sort5.30ident_100slen.d"
node_file <- "rp_nodes.xlsx"
output_csv <- "subsampling_modularity_results.csv"
output_plot <- "modularity_robustness_comparison.pdf"

n_replicates <- 100
sampling_size <- 336 # Size of the smallest major phylum

# 3. LOAD AND PREPARE ORIGINAL DATA ----
# --------------------------------------
print("Loading original node and edge data...")

# Load edges (assuming tab-delimited, no header)
rp_edges <- read.delim(edge_file, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
colnames(rp_edges) <- c("from", "to")

# Load nodes
rp_nodes <- read_excel(node_file)
# Ensure the first column is named 'name' for igraph
colnames(rp_nodes)[1] <- "name"

# Check for Phylum column
phylum_column_name <- "Phylum"
if (!phylum_column_name %in% names(rp_nodes)) {
  stop(paste("Node file must contain a column named '", phylum_column_name, "' for stratification."))
}

# --- CRITICAL DATA CLEANING STEP ---
# Ensure all nodes in the edge list actually exist in the node list.
valid_node_ids <- rp_nodes$name
print(paste("Original number of edges:", nrow(rp_edges)))

rp_edges_clean <- rp_edges %>%
  filter(from %in% valid_node_ids & to %in% valid_node_ids)

print(paste("Number of edges after cleaning:", nrow(rp_edges_clean)))


# --- Calculate Modularity of the ORIGINAL FULL Network ---
print("Building and analyzing the full original network...")
original_graph <- graph_from_data_frame(rp_edges_clean, directed = FALSE, vertices = rp_nodes)
original_graph_simple <- simplify(original_graph, remove.multiple = TRUE, remove.loops = TRUE)

# Use Walktrap community detection
original_communities <- cluster_walktrap(original_graph_simple, weights = NA)
original_modularity <- modularity(original_graph_simple, membership(original_communities))

print(paste("Modularity of the full, original network:", round(original_modularity, 4)))


# 4. PERFORM THE SUBSAMPLING EXPERIMENT ----
# ------------------------------------------
print(paste("Starting subsampling experiment with", n_replicates, "replicates..."))

modularity_results <- numeric(n_replicates)

# Check node counts per phylum
phylum_counts <- rp_nodes %>% count(!!sym(phylum_column_name))
print("Node counts per phylum:")
print(phylum_counts)

for (i in 1:n_replicates) {
  
  # --- Step A: Stratified Subsampling ---
  subsampled_nodes_df <- rp_nodes %>%
    group_by(!!sym(phylum_column_name)) %>%
    sample_n(size = min(n(), sampling_size)) %>%
    ungroup()
  
  subsampled_node_ids <- subsampled_nodes_df$name
  
  # --- Step B: Create Sub-Network ---
  subsampled_edges_df <- rp_edges_clean %>%
    filter(from %in% subsampled_node_ids & to %in% subsampled_node_ids)
  
  # --- Step C: Build and Analyze ---
  sub_graph <- graph_from_data_frame(subsampled_edges_df, directed = FALSE, vertices = subsampled_nodes_df)
  sub_graph_simple <- simplify(sub_graph, remove.multiple = TRUE, remove.loops = TRUE)
  
  if (gorder(sub_graph_simple) > 0 && gsize(sub_graph_simple) > 0) {
    sub_communities <- cluster_walktrap(sub_graph_simple, weights = NA)
    modularity_results[i] <- modularity(sub_graph_simple, membership(sub_communities))
  } else {
    modularity_results[i] <- NA
  }

  if (i %% 10 == 0) print(paste("Completed replicate", i))
}


# 5. SUMMARIZE AND SAVE RESULTS ----
# ----------------------------------
print("Subsampling experiment complete.")

mean_modularity <- mean(modularity_results, na.rm = TRUE)
sd_modularity <- sd(modularity_results, na.rm = TRUE)

summary_df <- data.frame(
  Analysis = c("Original Full Network", "Subsampled Balanced Network"),
  Mean_Modularity = c(original_modularity, mean_modularity),
  SD_Modularity = c(NA, sd_modularity),
  Number_of_Replicates = c(NA, n_replicates)
)

print("--- FINAL RESULTS ---")
print(summary_df)

# Save raw data
results_df <- data.frame(replicate = 1:n_replicates, modularity = modularity_results)
write.csv(results_df, file = output_csv, row.names = FALSE)
print(paste("Raw results saved to:", output_csv))


# 6. VISUALIZE THE RESULTS ----
# ----------------------------
print("Generating visualization...")

# Create the plot
comparison_plot <- ggplot(results_df, aes(x = "Subsampled Balanced", y = modularity)) +
  # Raincloud elements
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.2, 
               point_colour = NA, fill = "lightblue", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_point(shape = 21, color = "black", fill = "gray", alpha = 0.3, 
             position = position_jitter(width = 0.05, seed = 123)) +
  
  # Original value reference
  geom_hline(yintercept = original_modularity, linetype = "dashed", color = "red", linewidth = 1) +
  geom_point(aes(x = "Original Full", y = original_modularity), 
             shape = 23, fill = "red", color = "black", size = 5) +
  
  # Labels
  annotate("text", x = "Original Full", y = original_modularity, 
           label = paste("Original =", round(original_modularity, 3)), 
           vjust = -1.5, color = "red", fontface = "bold") +
  annotate("text", x = "Subsampled Balanced", y = max(modularity_results, na.rm=TRUE), 
           label = paste("Mean =", round(mean_modularity, 3), "±", round(sd_modularity, 3)),
           vjust = -1.5, fontface = "bold") +

  # Styling
  scale_x_discrete(name = "Network Type") +
  scale_y_continuous(name = "Modularity Score", 
                     limits = c(min(original_modularity, min(modularity_results)) - 0.02, 
                                max(original_modularity, max(modularity_results)) + 0.02)) +
  ggtitle("Robustness of Network Modularity to Uneven Sampling") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Print and Save
print(comparison_plot)
ggsave(output_plot, plot = comparison_plot, width = 8, height = 6)
print(paste("Plot saved to:", output_plot))