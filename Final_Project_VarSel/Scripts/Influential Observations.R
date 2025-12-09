# Influential Observations

# ==========================================
# Influence Diagnostics for Final Models
# ==========================================

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("car"))       install.packages("car")

library(tidyverse)
library(car)

# ---------- Paths ----------
base_path   <- "~/Math 5387 Regression Analysis/Final Project Fall 2025"
data_path   <- file.path(base_path, "Data")
vs_path     <- file.path(base_path, "Results", "Variable_Selection")
infl_path   <- file.path(base_path, "Results", "Influence")

dir.create(infl_path, recursive = TRUE, showWarnings = FALSE)

# ---------- Load datasets ----------
datasets_files <- list.files(data_path, pattern = "\\.RData$", full.names = TRUE)

datasets <- lapply(datasets_files, function(f) {
  e <- new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]
})

names(datasets) <- datasets_files |>
  basename() |>
  tools::file_path_sans_ext() |>
  gsub("[^A-Za-z0-9]+", "_", .)

# ---------- Load selected models from variable selection ----------
selected_models_file <- file.path(vs_path, "selected_models.rds")
if (!file.exists(selected_models_file)) {
  stop("selected_models.rds not found. Make sure you saved selected_vars_list in the variable selection script.")
}
selected_models <- readRDS(selected_models_file)

# Expected methods (names in selected_models[[i]])
method_order <- c(
  "backward_p", "forward_p", "stepwise_p",
  "backward_AIC", "forward_AIC", "stepwise_AIC",
  "subsets_AIC", "subsets_BIC", "subsets_AdjR2",
  "subsets_RMSE", "subsets_Cp"
)

# ---------- Helper: build model formula for a given method ----------
get_model_formula <- function(method_result, response_var = "risk") {
  # If the method already has a full formula, use it
  if (!is.null(method_result$final_model)) {
    return(method_result$final_model)
  }
  
  # Otherwise, try to build from its variable list (subsets methods)
  if (!is.null(method_result$df) &&
      "variable" %in% names(method_result$df) &&
      nrow(method_result$df) > 0) {
    
    vars <- unique(method_result$df$variable)
    rhs  <- paste(vars, collapse = " + ")
    return(as.formula(paste(response_var, "~", rhs)))
  }
  
  # If nothing to build, return NULL
  return(NULL)
}

# ---------- Helper: run influence diagnostics for a fitted model ----------
run_influence_diagnostics <- function(fit, dataset_name, method_name, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Basic quantities
  n <- nobs(fit)
  p <- length(coef(fit))      # includes intercept
  
  hat_vals   <- hatvalues(fit)
  cooks_vals <- cooks.distance(fit)
  rstud      <- rstudent(fit)
  rstand     <- rstandard(fit)
  dfb        <- dfbetas(fit)
  
  # --- Leverage threshold ---
  # Typical rule: h_ii > 3p/n considered high leverage
  lev_thresh <- 3 * p / n
  high_lev_idx <- which(hat_vals > lev_thresh)
  
  lev_df <- tibble(
    obs   = seq_along(hat_vals),
    hat   = hat_vals,
    cooks = cooks_vals,
    rstud = rstud,
    rstand = rstand,
    high_leverage = hat_vals > lev_thresh
  )
  
  write.csv(lev_df, file.path(out_dir, "influence_measures.csv"),
            row.names = FALSE)
  
  # Also save just flagged points for quick scanning
  flagged_df <- lev_df |>
    filter(high_leverage | abs(rstud) > 2 | cooks > (4 / n))
  
  write.csv(flagged_df, file.path(out_dir, "flagged_observations.csv"),
            row.names = FALSE)
  
  # --- Bonferroni outlier test on studentized residuals ---
  out_bonf <- tryCatch(
    outlierTest(fit),
    error = function(e) NULL
  )
  
  if (!is.null(out_bonf)) {
    bonf_file <- file.path(out_dir, "bonferroni_outliers.csv")
    # outlierTest returns an object with rows = outliers
    bonf_df <- as.data.frame(out_bonf)
    bonf_df$obs <- as.numeric(rownames(bonf_df))
    rownames(bonf_df) <- NULL
    write.csv(bonf_df, bonf_file, row.names = FALSE)
  }
  
  # ----------------------------
  # Plots
  # ----------------------------
  
  # 1) Half-normal & index plots of hat values
  png(file.path(out_dir, "hat_halfnormal.png"), width = 800, height = 600)
  qqPlot(abs(hat_vals),
         main = paste("Half-normal Plot of |hat| –", dataset_name, method_name),
         ylab = "|hat values|", id = FALSE)
  dev.off()
  
  png(file.path(out_dir, "hat_index.png"), width = 800, height = 600)
  plot(hat_vals, type = "h",
       main = paste("Index Plot of hat values –", dataset_name, method_name),
       xlab = "Observation index", ylab = "hat")
  abline(h = lev_thresh, col = "red", lty = 2)
  dev.off()
  
  # 2) Half-normal & index plots of Cook's distance
  png(file.path(out_dir, "cooks_halfnormal.png"), width = 800, height = 600)
  qqPlot(abs(cooks_vals),
         main = paste("Half-normal Plot of |Cook's D| –", dataset_name, method_name),
         ylab = "|Cook's distance|", id = FALSE)
  dev.off()
  
  png(file.path(out_dir, "cooks_index.png"), width = 800, height = 600)
  plot(cooks_vals, type = "h",
       main = paste("Index Plot of Cook's distance –", dataset_name, method_name),
       xlab = "Observation index", ylab = "Cook's D")
  abline(h = 4 / n, col = "red", lty = 2)
  dev.off()
  
  # 3) DFBETAs plots for each predictor
  png(file.path(out_dir, "dfbetas_plots.png"), width = 900, height = 900)
  dfbetaPlots(fit,
              main = paste("DFBETAs –", dataset_name, method_name))
  dev.off()
  
  # 4) Influence plot (std residuals vs leverage with Cook's D bubbles)
  png(file.path(out_dir, "influence_plot.png"), width = 800, height = 600)
  influencePlot(fit,
                main = paste("Influence Plot –", dataset_name, method_name),
                sub = "", id.method = "identify", 
                id.n = 0)  # don't auto-label; keeps plot clean
  dev.off()
  
  # 5) Standardized residuals vs leverage (with Cook's D contours)
  png(file.path(out_dir, "stdres_vs_leverage.png"), width = 800, height = 600)
  plot(hat_vals, rstand,
       xlab = "Leverage (hat)", ylab = "Standardized Residuals",
       main = paste("Std Residuals vs Leverage –", dataset_name, method_name))
  abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "black", "red"))
  
  # Overlay Cook's distance contours
  # Cook's D approx lines: D = const
  # D = (r^2 * h) / [p * (1 - h)]
  # We'll just identify big values via point size:
  points(hat_vals, rstand,
         cex = pmin(1 + 10 * cooks_vals, 5),
         col = rgb(0, 0, 1, 0.5))
  dev.off()
  
  # Return summary metrics for dashboard
  return(list(
    n_high_leverage = sum(high_lev_idx),
    n_big_cooksD    = sum(cooks_vals > (4/n)),
    n_big_rstudent  = sum(abs(rstud) > 2),
    n_bonf          = if (!is.null(out_bonf)) nrow(as.data.frame(out_bonf)) else 0
  ))
  
}

influence_path <- file.path(base_path, "Results", "Influential_Observations")
dir.create(influence_path, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(datasets)) {
  
  df_name <- names(datasets)[i]
  df      <- datasets[[i]]
  
  cat("\n=== Influential Observation Diagnostics for Dataset:", df_name, "===\n")
  
  # Uniqueness table for this dataset
  uniq <- model_uniqueness[[df_name]]
  canonical_map <- uniq$canonical      # method → canonical method
  keys          <- uniq$model_keys     # method → collapsed formula string
  
  dataset_dir <- file.path(influence_path, df_name)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (method in method_order) {
    
    # Skip methods that didn’t produce a model
    model_res <- selected_models[[df_name]][[method]]
    if (is.null(model_res)) next
    
    # Find the canonical representative for this method
    canonical_method <- canonical_map[[method]]
    
    # Only run diagnostics ONCE for each unique model
    # i.e., ONLY when this method IS the canonical representative
    if (canonical_method != method) {
      cat("  -> Skipping ", method,
          " (duplicate of ", canonical_method, ")\n", sep = "")
      next
    }
    
    # Build model formula
    model_formula <- get_model_formula(model_res, response_var = "risk")
    if (is.null(model_formula)) {
      warning("Could not extract formula for ", df_name, " method ", method)
      next
    }
    
    # Fit the model
    fit <- tryCatch(
      lm(model_formula, data = df),
      error = function(e) {
        warning("Model failed for dataset ", df_name, " method ", method,
                ": ", e$message)
        NULL
      }
    )
    if (is.null(fit)) next
    
    # Output directory for this (unique) model
    method_dir <- file.path(dataset_dir, canonical_method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat("  -> Running influential observation diagnostics for method: ",
        canonical_method, "\n")
    
    run_influence_diagnostics(
      fit          = fit,
      data         = df,
      dataset_name = df_name,
      method_name  = canonical_method,
      out_dir      = method_dir
    )
  }
}

# Convert dashboard_list to a clean data.frame
dashboard_df <- bind_rows(dashboard_list, .id = "Method")

dashboard_file <- file.path(ds_infl_dir, "influence_dashboard.csv")
write.csv(dashboard_df, dashboard_file, row.names = FALSE)

cat("  Dashboard saved for dataset:", ds_name, "\n")

