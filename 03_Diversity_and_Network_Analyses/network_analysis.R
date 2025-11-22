# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: network_analysis.R
#
# PURPOSE:
# This script constructs a Sequence Similarity Network (SSN) for viral RdRPs,
# visualizes it, and calculates its global modularity score.
#
# The network nodes represent viral clusters (e.g., at 50% identity), and
# edges represent significant BLASTp similarities between them.
#
# INPUT:
#   - Edge list: A tab-delimited file (Source, Target) from BLASTp results.
#   - Node metadata: An Excel/CSV file with node attributes (e.g., size, host).
#
# OUTPUT:
#   - A network visualization plot.
#   - The global modularity score of the network.
#
# USAGE:
# Run this script in RStudio or via command line: Rscript 04_network_analysis_and_modularity.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. LOAD LIBRARIES ----
# ----------------------
if (!require("igraph")) install.packages("igraph")
if (!require("qgraph")) install.packages("qgraph")
if (!require("readxl")) install.packages("readxl")

library(igraph)
library(qgraph)
library(readxl)

# 2. SETUP PATHS ----
# -------------------
# !!! UPDATE THESE PATHS !!!
input_dir <- "path/to/network_data"
edge_file <- file.path(input_dir, "rp_edges_30ident_100slen.txt")
node_file <- file.path(input_dir, "rp_nodes.xlsx")

# 3. LOAD DATA ----
# -----------------
print("Loading network data...")

# Load Edge List
# Expected format: Column 1 = Source ID, Column 2 = Target ID
rp_edges <- read.delim(edge_file, header = FALSE)

# Load Node Metadata
# Expected format: Column 1 = ID, followed by attributes (e.g., 'clstr50_num90', 'color_host')
rp_nodes <- read_excel(node_file)


# 4. CONSTRUCT AND SIMPLIFY NETWORK ----
# --------------------------------------
print("Constructing graph...")

# Create graph object
graph <- graph_from_data_frame(rp_edges, directed = FALSE, vertices = rp_nodes)

# Simplify graph (remove self-loops and multiple edges)
graph_simple <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

print(paste("Graph constructed with", vcount(graph_simple), "nodes and", ecount(graph_simple), "edges."))


# 5. CALCULATE LAYOUT ----
# ------------------------
print("Calculating Fruchterman-Reingold layout...")
# Using qgraph's layout for better handling of large networks
e <- get.edgelist(graph_simple, names = FALSE)
l <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(graph_simple),
                                       area = 8 * (vcount(graph_simple)^2),
                                       repulse.rad = (vcount(graph_simple)^3.1))


# 6. VISUALIZATION PARAMETERS ----
# --------------------------------
print("Setting visual attributes...")

# Node attributes
# Set node size based on the number of sequences in the cluster (log scale)
# 'clstr50_num90' must be a column in your node metadata
V(graph_simple)$size <- log(0.1 + V(graph_simple)$clstr50_num90) / 1.5
V(graph_simple)$frame.color <- NA
V(graph_simple)$label <- NA  # Remove labels for global view
V(graph_simple)$label.cex <- 0.01

# Color nodes by host group (assuming 'color_host' column in metadata)
V(graph_simple)$color <- V(graph_simple)$color_host

# Edge attributes
E(graph_simple)$color <- "lightgrey"
E(graph_simple)$width <- 0.01


# 7. PLOT NETWORK ----
# --------------------
print("Plotting network...")
# Save to PDF
pdf("network_plot.pdf", width = 10, height = 10)
plot(graph_simple, layout = l)
dev.off()


# 8. COMMUNITY DETECTION AND MODULARITY ----
# ------------------------------------------
print("Calculating modularity...")

# Detect communities using the Walktrap algorithm
graph_cwalk <- cluster_walktrap(graph_simple, weights = NA, steps = 12)

# Calculate global modularity score
mod_score <- modularity(graph_simple, membership(graph_cwalk))

print(paste("Global Modularity Score:", round(mod_score, 4)))

# Optional: Visualize communities
# plot(graph_cwalk, graph_simple, layout = l)