# Error Assumptions
# ==========================================
# Error Assumptions Diagnostics for Final Models
# ==========================================

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("car"))       install.packages("car")

library(tidyverse)
library(car)

# ---------- Paths ----------
base_path   <- "~/Math 5387 Regression Analysis/Final Project Fall 2025"
data_path   <- file.path(base_path, "Data")
vs_path     <- file.path(base_path, "Results", "Variable_Selection")
err_path    <- file.path(base_path, "Results", "Error_Assumptions")

dir.create(err_path, recursive = TRUE, showWarnings = FALSE)

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

# Methods to iterate over (same order as before)
method_order <- c(
  "backward_p", "forward_p", "stepwise_p",
  "backward_AIC", "forward_AIC", "stepwise_AIC",
  "subsets_AIC", "subsets_BIC", "subsets_AdjR2",
  "subsets_RMSE", "subsets_Cp"
)

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

error_assumptions_diagnostics <- function(fit, data, dataset_name, method_name, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  # Residuals & fitted
  resid_used   <- residuals(fit)
  std_resid    <- rstandard(fit)
  fitted_vals  <- fitted(fit)
  n            <- length(resid_used)
  
  # For variance ratio tests
  mf          <- model.frame(fit)               # rows used in fit
  term_labels <- attr(terms(fit), "term.labels")
  
  # ================================
  # 1) Mean zero – Pearson vs fitted
  # ================================
  pearson_resid <- residuals(fit, type = "pearson")
  
  png(file.path(out_dir, "pearson_vs_fitted.png"), width = 800, height = 600)
  plot(fitted_vals, pearson_resid,
       xlab = "Fitted values",
       ylab = "Pearson residuals",
       main = paste("Pearson Residuals vs Fitted –", dataset_name, method_name))
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  
  # ================================
  # 2) Constant variance
  #    - Scale–location plot
  #    - Variance ratio F-tests
  # ================================
  
  # Scale-location: sqrt(|std residuals|) vs fitted
  png(file.path(out_dir, "scale_location.png"), width = 800, height = 600)
  plot(fitted_vals, sqrt(abs(std_resid)),
       xlab = "Fitted values",
       ylab = expression(sqrt("|Standardized residuals|")),
       main = paste("Scale-Location Plot –", dataset_name, method_name))
  dev.off()
  
  # Variance ratio F-tests
  numeric_terms <- term_labels[!grepl(":", term_labels)]
  numeric_terms <- numeric_terms[numeric_terms %in% names(mf)]
  numeric_terms <- numeric_terms[sapply(mf[numeric_terms], is.numeric)]
  
  vr_df <- data.frame()
  
  for (v in numeric_terms) {
    x   <- mf[[v]]
    mid <- median(x, na.rm = TRUE)   # “midpoint”
    
    g1 <- resid_used[x <= mid]
    g2 <- resid_used[x > mid]
    
    if (length(g1) > 1 && length(g2) > 1) {
      vt <- var.test(g1, g2)
      vr_df <- rbind(
        vr_df,
        data.frame(
          method   = method_name,
          variable = v,
          midpoint = mid,
          F_stat   = as.numeric(vt$statistic),
          df1      = as.numeric(vt$parameter[1]),
          df2      = as.numeric(vt$parameter[2]),
          p_value  = vt$p.value
        )
      )
    }
  }
  results$variance_ratio <- vr_df
  
  # ================================
  # 3) Normality
  #    - QQ plot
  #    - Shapiro–Wilk
  # ================================
  
  png(file.path(out_dir, "qqplot_residuals.png"), width = 800, height = 600)
  qqPlot(resid_used,
         main = paste("QQ Plot of Residuals –", dataset_name, method_name),
         ylab = "Residuals", id = FALSE)
  dev.off()
  
  sh <- shapiro.test(resid_used)
  sh_df <- data.frame(
    method  = method_name,
    W       = as.numeric(sh$statistic),
    p_value = sh$p.value
  )
  results$shapiro <- sh_df
  
  # ================================
  # 4) Correlated errors
  #    - Durbin–Watson
  #    - Serial lag-1 plot + correlation
  # ================================
  
  dw <- durbinWatsonTest(fit)
  dw_df <- data.frame(
    method  = method_name,
    DW      = as.numeric(dw$dw),
    p_value = dw$p,
    lag     = dw$lags[1]
  )
  results$durbin_watson <- dw_df
  
  # Serial correlation plot: e_(i+1) vs e_i
  if (length(resid_used) > 1) {
    e_i   <- head(resid_used, n - 1)
    e_ip1 <- tail(resid_used, n - 1)
    lag1  <- cor(e_i, e_ip1)
    
    png(file.path(out_dir, "serial_correlation_ei_ei1.png"),
        width = 800, height = 600)
    plot(e_i, e_ip1,
         xlab = expression(hat(epsilon)[i]),
         ylab = expression(hat(epsilon)[i+1]),
         main = paste("Lag-1 Residual Plot –", dataset_name, method_name))
    abline(h = 0, v = 0, col = grey(0.75))
    dev.off()
  } else {
    lag1 <- NA
  }
  
  lag_df <- data.frame(
    method = method_name,
    lag1_correlation = lag1
  )
  results$lag1_corr <- lag_df
  
  return(results)
}


error_path <- file.path(base_path, "Results", "Error_Assumptions")
dir.create(error_path, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(datasets)) {
  
  ds_name <- names(datasets)[i]
  df      <- datasets[[i]]
  
  cat("\n=== Error Assumption Diagnostics for Dataset:", ds_name, "===\n")
  
  # Uniqueness mapping for this dataset
  uniq <- model_uniqueness[[ds_name]]
  canonical_map <- uniq$canonical      # method → canonical representative
  keys          <- uniq$model_keys     # method → compact model string
  
  # Dataset-level output folder
  ds_dir <- file.path(error_path, ds_name)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Dataset-level summary storage
  vr_all  <- data.frame()   # variance ratio F-tests
  sh_all  <- data.frame()   # Shapiro-Wilk
  dw_all  <- data.frame()   # Durbin-Watson
  lag_all <- data.frame()   # lag-1 correlations
  
  for (method in method_order) {
    
    model_res <- selected_models[[ds_name]][[method]]
    if (is.null(model_res)) next
    
    canonical_method <- canonical_map[[method]]
    
    # Only evaluate model for its first (canonical) occurrence
    if (canonical_method != method) {
      cat("  -> Skipping ", method,
          " (duplicate of ", canonical_method, ")\n", sep = "")
      next
    }
    
    # Extract formula
    form <- get_model_formula(model_res, response_var = "risk")
    if (is.null(form)) {
      warning("Could not extract formula for dataset ", ds_name, ", method ", method)
      next
    }
    
    # Fit model
    fit <- tryCatch(
      lm(form, data = df),
      error = function(e) {
        warning("Model failed for dataset ", ds_name, " method ", method, 
                ": ", e$message)
        NULL
      }
    )
    if (is.null(fit)) next
    
    cat("  -> Running diagnostics for unique model: ", canonical_method, "\n")
    
    # Folder for plots for this unique model
    method_dir <- file.path(ds_dir, canonical_method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Plot diagnostics + return summary tables
    diag_res <- error_assumptions_diagnostics(
      fit          = fit,
      data         = df,
      dataset_name = ds_name,
      method_name  = canonical_method,
      out_dir      = method_dir
    )
    
    # Append results to dataset-level summaries
    vr_all  <- rbind(vr_all,  diag_res$variance_ratio)
    sh_all  <- rbind(sh_all,  diag_res$shapiro)
    dw_all  <- rbind(dw_all,  diag_res$durbin_watson)
    lag_all <- rbind(lag_all, diag_res$lag1_corr)
  }
  
  # ------------------------------------
  # Write a single combined CSV per dataset
  # ------------------------------------
  
  summary_file <- file.path(ds_dir, paste0(ds_name, "_error_assumptions_summary.csv"))
  con <- file(summary_file, "w")
  
  writeLines("Variance Ratio F-tests", con)
  if (nrow(vr_all) > 0) {
    write.table(vr_all, con, sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
  } else {
    writeLines("None\n", con)
  }
  
  writeLines("\n\nShapiro-Wilk Normality Tests", con)
  if (nrow(sh_all) > 0) {
    write.table(sh_all, con, sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
  } else {
    writeLines("None\n", con)
  }
  
  writeLines("\n\nDurbin-Watson Tests", con)
  if (nrow(dw_all) > 0) {
    write.table(dw_all, con, sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
  } else {
    writeLines("None\n", con)
  }
  
  writeLines("\n\nLag-1 Residual Correlations", con)
  if (nrow(lag_all) > 0) {
    write.table(lag_all, con, sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
  } else {
    writeLines("None\n", con)
  }
  
  close(con)
  
  cat("  ✔ Summary saved:", summary_file, "\n")
}
