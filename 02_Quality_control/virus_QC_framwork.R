# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Script to Create Host Veracity QC Plots (Strip & Box)
#
# This script reads a single QC metrics file and generates the final
# figure comparing 'Trusted' vs. 'Other' data for only the two
# statistically significant metrics: contigs ratio and virus abundance.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Load libraries
# -------------------
library(tidyverse)
library(patchwork)
library(ggtext)
library(Cairo)

# 2. Setup paths, colors, and labels
# ------------------------------------
# !!! PLEASE UPDATE THIS PATH !!!
input_file <- file.path(base_path, "virus_QC_metrix.txt")
output_file <- file.path(base_path, "virus_QC_Figure.pdf")

color_palette_box <- c("Trusted" = "#ff4500", "Other" = "#808080")
color_palette_points <- c("Trusted" = "#ff4500", "Other" = "grey50")
group_order <- c("Trusted", "Other")


# 3. Load and prepare data
# --------------------------
print(paste("Loading and preparing data from:", input_file))
qc_data <- read_tsv(input_file, show_col_types = FALSE) %>%
  rename(
    log_rpm = log_RPM,
    contigs_ratio = contigs_Landplants_Percent,
    group = Trusted_data
  ) %>%
  mutate(
    group_label = factor(if_else(group == "T", "Trusted", "Other"), levels = group_order)
  ) %>%
  # Only select the columns we actually need
  select(Reference_virus_ID, log_rpm, contigs_ratio, group_label) %>%
  drop_na()

print("Data loaded and prepared successfully.")


# 4. Helper functions
# ---------------------
# Function for P-value annotation
get_p_value_label <- function(data, value_col) {
  test_result <- wilcox.test(as.formula(paste(value_col, "~ group_label")), 
                             data = data, alternative = "greater")
  p_value <- test_result$p.value
  if (p_value < 0.001) {
    "<i>P</i> < 0.001"
  } else {
    paste0("<i>P</i> = ", format(p_value, nsmall = 3, digits = 3))
  }
}

# Master function to create a column of plots
create_plot_column <- function(data, x_var, x_lab, plot_title, x_lim) {
  
  trusted_group_data <- data %>% filter(group_label == "Trusted")
  median_trusted <- median(trusted_group_data[[x_var]], na.rm = TRUE)
  quantiles_trusted <- quantile(trusted_group_data[[x_var]], probs = c(0.025, 0.975), na.rm = TRUE)
  
  annotation_text <- glue::glue(
    "Median: {format(median_trusted, nsmall=2, digits=2)}<br>95% Range: [{format(quantiles_trusted[1], nsmall=2, digits=2)}, {format(quantiles_trusted[2], nsmall=2, digits=2)}]"
  )
  
  # --- Strip Plot ---
  p_strip <- ggplot(data, aes(x = .data[[x_var]], y = 0)) +
    geom_jitter(data = . %>% filter(group_label == "Other"), 
                color = color_palette_points["Other"], alpha = 0.4, height = 0.35, size = 1, shape = 16) +
    geom_jitter(data = . %>% filter(group_label == "Trusted"),
                color = color_palette_points["Trusted"], alpha = 0.9, height = 0.35, size = 1, shape = 16) +
    annotate("rect", xmin = quantiles_trusted[1], xmax = quantiles_trusted[2], ymin = -Inf, ymax = Inf, fill = "#ff9999", alpha = 0.15) +
    geom_vline(xintercept = median_trusted, color = "#ff4500", linetype = "dashed", alpha = 0.7) +
    annotate("richtext", x = -Inf, y = Inf, label = annotation_text,
             hjust = -0.05, vjust = 1.05, size = 4) + # Slightly larger font
    scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0.05, 0.05))) +
    labs(title = plot_title, x = NULL, y = NULL) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), # Slightly larger
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_line(color = "grey92")
    )

  # --- Box Plot ---
  p_box <- ggplot(data, aes(x = .data[[x_var]], y = group_label, fill = group_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) + 
    annotate("richtext", x = -Inf, y = Inf, label = get_p_value_label(data, x_var),
             hjust = -0.1, vjust = 1.2, size = 4) + # Slightly larger font
    scale_fill_manual(values = color_palette_box) +
    scale_y_discrete(limits = rev(levels(data$group_label))) +
    scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = x_lab, y = NULL) +
    theme_bw(base_family = "Arial") +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "bold"), # Slightly larger
      axis.title.x = element_text(size = 14), # Slightly larger
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  plot_column <- p_strip / p_box + plot_layout(heights = c(1.25, 0.75))
  
  return(plot_column)
}


# 5. Generate and assemble final figure
# ---------------------------------------
print("Generating plot columns for the two validated metrics...")
col_contigs <- create_plot_column(qc_data, "contigs_ratio", "contigs ratio (%)", "Contigs plants/eukaryotes ratio", c(0, 100))
col_abundance <- create_plot_column(qc_data, "log_rpm", "log10(RPM)", "virus abundance", c(-1.5, 6))

# Assemble the final 2x2 figure
final_plot <- col_contigs | col_abundance

# Add the final A, B, C, D panel labels
final_plot_labeled <- final_plot + 
  plot_annotation(tag_levels = list(c('A', 'C', 'B', 'D'))) &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Arial"))


# 6. Save the plot
# ------------------
print("Saving final figure...")
ggsave(
  output_file,
  plot = final_plot_labeled,
  width = 12,
  height = 7,
  units = "in",
  device = cairo_pdf,
  dpi = 600
)
print(paste("Script finished. Final QC figure saved to:", output_file))