# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Script for Rarefaction Analysis based on Library Count (Original Logic)
# Version 8.0
#
# This script performs rarefaction analyses using the number of libraries
# as the sampling unit. It reads long-format discovery data, converts it
# to the necessary incidence matrices, and runs iNEXT.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. SETUP: LOAD LIBRARIES AND DEFINE PARAMETERS ----
# --------------------------------------------------
library(tidyverse)
library(iNEXT)
library(ggplot2)

# --- Parameters ---
N_BOOT <- 200

# --- Input File Paths ---
PREPROCESSED_DIR <- "preprocessed_for_inext_long" 

# --- Analysis Configurations ---
config_E <- list(
  name = "Analysis_E_by_AAI",
  files = c("AAI_0.3" = "AAI_0.3.long.csv", "AAI_0.5" = "AAI_0.5.long.csv",
            "AAI_0.7" = "AAI_0.7.long.csv", "AAI_0.9" = "AAI_0.9.long.csv")
)

config_F <- list(
  name = "Analysis_F_by_Plant_Group",
  files = c("Bryophytes" = "Bryophytes.long.csv", "Lycophytes" = "Lycophytes.long.csv",
            "Ferns" = "Ferns.long.csv", "Gymnosperms" = "Gymnosperms.long.csv",
            "Angiosperms" = "Angiosperms.long.csv")
)

config_G <- list(
  name = "Analysis_G_by_Virus_Phylum",
  files = c("Duplornaviricota" = "Duplornaviricota.long.csv",
            "Kitrinoviricota" = "Kitrinoviricota.long.csv",
            "Lenarviricota" = "Lenarviricota.long.csv",
            "Negarnaviricota" = "Negarnaviricota.long.csv",
            "Pisuviricota" = "Pisuviricota.long.csv")
)

all_analyses <- list(config_E, config_F, config_G)


# 2. HELPER FUNCTION TO CREATE INCIDENCE MATRIX ----
# --------------------------------------------------
# This function takes long data (library_id, virus_cluster) and converts it
# to a wide incidence matrix (rows=viruses, cols=libraries).
# Crucially, it treats each library ID as a single sampling unit.
create_incidence_matrix <- function(long_df) {
  long_df %>%
    distinct(library_id, virus_cluster) %>% # Ensure unique hits
    mutate(value = 1) %>%
    pivot_wider(names_from = library_id, values_from = value, values_fill = 0) %>%
    column_to_rownames(var = "virus_cluster")
}


# 3. CORE ANALYSIS LOOP ----
# --------------------------
run_rarefaction_analysis <- function(config) {
  
  analysis_name <- config$name
  group_files <- config$files
  
  print(paste("--- Starting:", analysis_name, "---"))
  
  inext_input_list <- list()
  
  for (group_name in names(group_files)) {
    file_name <- group_files[[group_name]]
    file_path <- file.path(PREPROCESSED_DIR, file_name)
    
    if (file.exists(file_path)) {
      print(paste("  Loading:", group_name))
      
      # Read the long-format file
      # It has columns: virus_cluster, library_id (and possibly others from Python script)
      # We only need those two.
      discoveries_long <- read_csv(file_path, show_col_types = FALSE)
      
      # Create the standard incidence matrix
      # Dimensions: [Number of Viruses] x [Number of Libraries]
      if (nrow(discoveries_long) > 0) {
        wide_matrix <- create_incidence_matrix(discoveries_long)
        inext_input_list[[group_name]] <- wide_matrix
      } else {
        print(paste("    Warning: No data found in", file_name))
      }
    } else {
      print(paste("    Error: File not found:", file_path))
    }
  }
  
  # --- Run iNEXT ---
  if (length(inext_input_list) > 0) {
    print("  Running iNEXT calculation...")
    # datatype="incidence_raw" tells iNEXT that columns are sampling units (libraries)
    output_inext <- iNEXT(inext_input_list, datatype = "incidence_raw", nboot = N_BOOT)
    
    # --- Create Plot ---
    # type=1 is the sample-size-based rarefaction curve
    plot_final <- ggiNEXT(output_inext, type = 1) + 
      labs(title = analysis_name,
           x = "Sampling Effort (Number of Transcriptomes)",
           y = "Number of Viral Clusters") +
      theme_bw(base_size = 14)
    
    return(list(plot = plot_final, output = output_inext))
  } else {
    return(NULL)
  }
}


# 4. EXECUTE ALL ANALYSES ----
# ----------------------------
for (config in all_analyses) {
  result <- run_rarefaction_analysis(config)
  
  if (!is.null(result)) {
    final_plot <- result$plot
    print(final_plot)
    
    # Save the main plot
    ggsave(paste0(config$name, "_rarefaction_original_metric.pdf"), final_plot, width = 10, height = 8)
    
    # Create Zoom-in plot for Analysis F (Plant Groups)
    if (config$name == "Analysis_F_by_Plant_Group") {
      plot_zoom <- final_plot +
        coord_cartesian(xlim = c(0, 2000), ylim = c(0, 400)) + # Adjust ylim based on your data range
        labs(title = paste(config$name, "(Zoomed-in View)"))
      
      print(plot_zoom)
      ggsave(paste0(config$name, "_rarefaction_original_metric_zoom.pdf"), plot_zoom, width = 10, height = 8)
    }
  }
}

print("--- All analyses complete. ---")