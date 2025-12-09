visualize_dataset <- function(data, remaining_vars, dataset_name, output_path) {
  
  # Folder
  plot_dir <- file.path(output_path, dataset_name)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Helpers
  clean_label <- function(x) gsub("_", " ", tools::toTitleCase(x))
  clean_filename <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  
  # Split variables into types
  numeric_vars  <- remaining_vars[sapply(data[remaining_vars], is.numeric)]
  factor_vars   <- remaining_vars[sapply(data[remaining_vars], is.factor)]
  
  # --------------------------------------------
  # 1️⃣ FACTOR VARIABLE VISUALIZATIONS
  # --------------------------------------------
  for (var in factor_vars) {
    
    ## Barplot
    bar_file <- file.path(plot_dir, paste0(clean_filename(var), "_barplot.png"))
    g1 <- ggplot(data, aes(x = .data[[var]])) +
      geom_bar(fill = "#59A14F", color = "black") +
      labs(title = paste("Counts of", clean_label(var)),
           x = clean_label(var), y = "Count") +
      theme_minimal(base_size = 14)
    ggsave(bar_file, g1, width = 6, height = 4, dpi = 300)
    
    ## Boxplot vs Risk
    box_file <- file.path(plot_dir, paste0(clean_filename(var), "_vs_risk.png"))
    g2 <- ggplot(data, aes(x = .data[[var]], y = risk)) +
      geom_boxplot(fill = "#E15759", alpha = 0.7) +
      labs(title = paste("Risk by", clean_label(var)),
           x = clean_label(var), y = "Risk") +
      theme_minimal(base_size = 14)
    ggsave(box_file, g2, width = 6, height = 4, dpi = 300)
  }
  
  # --------------------------------------------
  # 2️⃣ NUMERIC VARIABLE VISUALIZATIONS
  # --------------------------------------------
  for (var in numeric_vars) {
    
    ## Univariate plot (density vs histogram)
    if (any(data[[var]] %% 1 != 0, na.rm = TRUE)) {
      # Float → density
      uni_file <- file.path(plot_dir, paste0(clean_filename(var), "_density.png"))
      g_uni <- ggplot(data, aes(x = .data[[var]])) +
        geom_density(fill = "#5DA5DA", alpha = 0.6, color = "black") +
        labs(title = paste("Density of", clean_label(var)),
             x = clean_label(var), y = "Density") +
        theme_minimal(base_size = 14)
    } else {
      # Integer → histogram
      uni_file <- file.path(plot_dir, paste0(clean_filename(var), "_hist.png"))
      g_uni <- ggplot(data, aes(x = .data[[var]])) +
        geom_histogram(bins = 30, fill = "#F28E2B", color = "black") +
        labs(title = paste("Histogram of", clean_label(var)),
             x = clean_label(var), y = "Count") +
        theme_minimal(base_size = 14)
    }
    ggsave(uni_file, g_uni, width = 6, height = 4, dpi = 300)
    
    ## Scatter vs Risk
    scatter_file <- file.path(plot_dir, paste0(clean_filename(var), "_vs_risk.png"))
    g_scatter <- ggplot(data, aes(x = .data[[var]], y = risk)) +
      geom_point(alpha = 0.5, color = "#4E79A7") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = paste(clean_label(var), "vs Risk"),
           x = clean_label(var), y = "Risk") +
      theme_minimal(base_size = 14)
    ggsave(scatter_file, g_scatter, width = 6, height = 4, dpi = 300)
  }
  return(TRUE)
}

# ============================================================
# Run Visualizations for Each Dataset After Collinearity Step
# ============================================================

library(tidyverse)

# ----- Paths -----
base_path <- "/cloud/project/Fall 2025/Final Project Fall 2025"
data_path <- file.path(base_path, "Data")
collinearity_path <- file.path(base_path, "Results", "Collinearity")
viz_output_path <- file.path(base_path, "Results", "Data_Visualization")
dir.create(viz_output_path, recursive = TRUE, showWarnings = FALSE)

# ----- File with remaining variables for ALL datasets -----
remaining_file <- file.path(collinearity_path, "Remaining_Variables_All_Datasets.txt")

if (!file.exists(remaining_file)) {
  stop("ERROR: Could not find 'remaining_variables_all_datasets.txt' in Collinearity folder.")
}

# ---- Read & Parse Remaining Variables File ----
raw <- readLines(remaining_file)

remaining_list <- list()
current_dataset <- NULL
collecting <- FALSE

for (line in raw) {
  
  stripped <- trimws(line)
  
  # Detect dataset header
  if (grepl("^Dataset:", stripped)) {
    # Extract dataset name
    current_dataset <- gsub("^Dataset:\\s*", "", stripped)
    current_dataset <- gsub("\\s+", "_", current_dataset)  # clean name (remove spaces)
    remaining_list[[current_dataset]] <- character()
    collecting <- FALSE
    next
  }
  
  # Begin collecting after "Remaining Variables (sorted):"
  if (grepl("^Remaining Variables", stripped)) {
    collecting <- TRUE
    next
  }
  
  # Collect variable lines like " - varname"
  if (collecting && grepl("^- ", stripped)) {
    var <- gsub("^- ", "", stripped)
    remaining_list[[current_dataset]] <- c(remaining_list[[current_dataset]], var)
    next
  }
  
  # Stop collecting when next ==== block appears
  if (grepl("^====", stripped)) {
    collecting <- FALSE
    next
  }
}

# ---- Load datasets from Data folder ----
datasets <- lapply(list.files(data_path, pattern = "\\.RData$", full.names = TRUE), function(f) {
  e <- new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]  # single object
})

dataset_names <- gsub("\\.RData$", "", list.files(data_path, pattern = "\\.RData$"))
dataset_names <- gsub("[^A-Za-z0-9]+", "_", dataset_names)
names(datasets) <- dataset_names

# ---- Visualization Output Path ----
vis_path <- file.path(dirname(data_path), "Results", "Data_Visualization")
dir.create(vis_path, recursive = TRUE, showWarnings = FALSE)

# ---- Run Visualization for Each Dataset ----
for (ds_name in names(datasets)) {
  
  if (!ds_name %in% names(remaining_list)) {
    message("WARNING: No remaining vars recorded for dataset: ", ds_name, ". Skipping.")
    next
  }
  
  remaining_vars <- remaining_list[[ds_name]]
  
  # Safety filtering so no missing variables cause failure
  remaining_vars <- remaining_vars[remaining_vars %in% names(datasets[[ds_name]])]
  
  if (length(remaining_vars) == 0) {
    message("WARNING: Dataset ", ds_name, " has zero remaining vars after filtering. Skipping.")
    next
  }
  
  message("Visualizing dataset: ", ds_name)
  
  visualize_dataset(
    data = datasets[[ds_name]],
    remaining_vars = remaining_vars,
    dataset_name = ds_name,
    output_path = vis_path
  )
}

message("✔ All visualizations completed and saved.")
