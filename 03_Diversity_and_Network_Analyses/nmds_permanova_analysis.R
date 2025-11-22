# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT: nmds_permanova_analysis.R
#
# PURPOSE:
# This script performs Non-metric Multidimensional Scaling (NMDS) analysis
# to visualize the dissimilarity in viral community composition among different
# eukaryote host groups. It also conducts Permutational Multivariate Analysis
# of Variance (PERMANOVA) to statistically test these differences.
#
# INPUT:
#   - otu.txt: A species abundance matrix (rows = samples, cols = viral OTUs)
#   - group.txt: A metadata file linking samples to host groups (cols = samples, sites)
#
# OUTPUT:
#   - An NMDS plot (PDF)
#   - PERMANOVA statistical results (CSV)
#
# USAGE:
# Run this script in RStudio or via command line: Rscript 08_nmds_permanova_analysis.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. LOAD LIBRARIES ----
# ----------------------
if (!require("vegan")) install.packages("vegan")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyverse")) install.packages("tidyverse")

library(vegan)
library(ggplot2)
library(tidyverse)

# 2. SETUP PATHS AND PARAMETERS ----
# ----------------------------------
# !!! UPDATE THESE PATHS !!!
input_dir <- "path/to/input_data"
output_dir <- "path/to/output_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

otu_file <- file.path(input_dir, "otu.txt")
group_file <- file.path(input_dir, "group.txt")

# 3. LOAD AND PREPARE DATA ----
# -----------------------------
print("Loading data...")

# Load OTU table (viral abundance)
# Assumes row names are samples
otu <- read.delim(otu_file, row.names = 1, header = TRUE, check.names = FALSE)

# Load Group metadata
group <- read.table(group_file, sep = '\t', header = TRUE)

# Ensure samples match between OTU and Group files
common_samples <- intersect(rownames(otu), group$samples)
otu <- otu[common_samples, ]
group <- group %>% filter(samples %in% common_samples)

# 4. NMDS ANALYSIS ----
# ---------------------
print("Calculating Bray-Curtis distances and running NMDS...")

# Calculate Bray-Curtis dissimilarity matrix
otu.distance <- vegdist(otu, method = 'bray')

# Run NMDS (k=2 dimensions)
# metaMDS automatically transforms data and finds the best solution
nmds_result <- metaMDS(otu.distance, k = 2, trymax = 100)

# Extract stress value
stress_value <- nmds_result$stress
print(paste("NMDS Stress:", round(stress_value, 4)))

# Prepare data for plotting
nmds_points <- as.data.frame(nmds_result$points)
nmds_points$samples <- rownames(nmds_points)
names(nmds_points)[1:2] <- c('NMDS1', 'NMDS2')

# Merge with group information
plot_data <- merge(nmds_points, group, by = "samples")

# 5. PLOT NMDS ----
# -----------------
print("Generating NMDS plot...")

# Define a custom color palette (optional)
my_colors <- c("#648c79", "#82aa72", "#d1e5aa", "#aa535e", 
               "#408fb4", "#f1c27a", "#e2844b", "#8182ae")

p_nmds <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = sites)) +
  theme_bw() +
  # Add points
  geom_point(size = 3, alpha = 0.8) +
  # Add confidence ellipses
  stat_ellipse(aes(fill = sites), geom = "polygon", level = 0.95, alpha = 0.2, show.legend = FALSE) +
  # Add center lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  # Custom colors
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  # Labels and Title
  labs(title = paste("NMDS of Viral Communities (Stress =", round(stress_value, 3), ")"),
       color = "Host Group") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Save plot
ggsave(file.path(output_dir, "NMDS_plot.pdf"), p_nmds, width = 8, height = 6)
print(p_nmds)


# 6. STATISTICAL ANALYSIS (PERMANOVA) ----
# ----------------------------------------
print("Running PERMANOVA tests...")

# --- A. Test for Homogeneity of Dispersions (BETADISPER) ---
# Ideally, groups should have similar dispersions for PERMANOVA to be strictly valid
dispersion <- betadisper(otu.distance, group$sites)
disp_test <- permutest(dispersion, permutations = 999)
print("Beta-dispersion test results:")
print(disp_test)

# --- B. Global PERMANOVA (ADONIS) ---
# Tests if there are ANY significant differences among groups
global_adonis <- adonis2(otu.distance ~ sites, data = group, permutations = 999, method = "bray")
print("Global PERMANOVA results:")
print(global_adonis)

# Save global results
write.csv(as.data.frame(global_adonis), file.path(output_dir, "permanova_global_results.csv"))


# --- C. Pairwise PERMANOVA ---
# Tests differences between specific pairs of host groups
print("Running Pairwise PERMANOVA...")

group_names <- unique(group$sites)
pairwise_results <- data.frame()

# Loop through all unique pairs
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    
    # Subset data for the current pair
    pair_groups <- c(as.character(group_names[i]), as.character(group_names[j]))
    subset_samples <- group %>% filter(sites %in% pair_groups)
    subset_otu <- otu[subset_samples$samples, ]
    
    # Run ADONIS on the subset
    # Note: Using adonis2 is generally preferred over adonis
    pair_adonis <- adonis2(subset_otu ~ sites, data = subset_samples, permutations = 999, method = 'bray')
    
    # Extract stats
    res_row <- data.frame(
      Group1 = group_names[i],
      Group2 = group_names[j],
      R2 = pair_adonis$R2[1],
      F_Model = pair_adonis$F[1],
      P_value = pair_adonis$`Pr(>F)`[1]
    )
    
    pairwise_results <- rbind(pairwise_results, res_row)
  }
}

# Adjust P-values for multiple comparisons (Bonferroni correction)
pairwise_results$P_adj <- p.adjust(pairwise_results$P_value, method = "bonferroni")

# Save pairwise results
print("Pairwise PERMANOVA results (Head):")
print(head(pairwise_results))
write.csv(pairwise_results, file.path(output_dir, "permanova_pairwise_results.csv"), row.names = FALSE)

print("--- Analysis Complete ---")