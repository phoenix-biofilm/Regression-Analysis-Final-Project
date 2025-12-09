# Prediction

# ============================================================
# Prediction Script (Run AFTER influential obs removal)
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
base_path <- "~/Math 5387 Regression Analysis/Final Project Fall 2025"
clean_data_path <- file.path(base_path, "Data_After_Influence_Removal")
vs_path   <- file.path(base_path, "Results", "Variable_Selection")
pred_path <- file.path(base_path, "Results", "Prediction")
dir.create(pred_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load datasets (already cleaned)
# ------------------------------------------------------------
dataset_files <- list.files(clean_data_path, pattern = "\\.RData$", full.names = TRUE)

datasets <- lapply(dataset_files, function(f) {
  e <- new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]
})

names(datasets) <- dataset_files |>
  basename() |>
  tools::file_path_sans_ext() |>
  gsub("[^A-Za-z0-9]+", "_", .)

# ------------------------------------------------------------
# Load variable selection results
# ------------------------------------------------------------
selected_models  <- readRDS(file.path(vs_path, "selected_models.rds"))
model_uniqueness <- readRDS(file.path(vs_path, "model_uniqueness_map.rds"))

method_order <- c(
  "backward_p", "forward_p", "stepwise_p",
  "backward_AIC", "forward_AIC", "stepwise_AIC",
  "subsets_AIC", "subsets_BIC", "subsets_AdjR2",
  "subsets_RMSE", "subsets_Cp"
)

# ------------------------------------------------------------
# Helper: extract model formula
# ------------------------------------------------------------
get_model_formula <- function(method_res, response_var = "risk") {
  if (!is.null(method_res$final_model)) return(method_res$final_model)
  
  if (!is.null(method_res$df) &&
      "variable" %in% names(method_res$df) &&
      nrow(method_res$df) > 0) {
    rhs <- paste(method_res$df$variable, collapse = " + ")
    return(as.formula(paste(response_var, "~", rhs)))
  }
  
  return(NULL)
}

# ------------------------------------------------------------
# Helper: extract numeric and factor variables used in model
# ------------------------------------------------------------
get_model_components <- function(fit, data) {
  terms_obj <- terms(fit)
  labels <- attr(terms_obj, "term.labels")
  
  main_terms <- labels[!grepl(":", labels)]
  main_terms <- main_terms[main_terms %in% names(data)]
  
  numeric_vars <- main_terms[sapply(data[main_terms], is.numeric)]
  factor_vars  <- main_terms[sapply(data[main_terms], is.factor)]
  
  # also detect interactions
  interactions <- labels[grepl(":", labels)]
  for (it in interactions) {
    parts <- strsplit(it, ":")[[1]]
    for (p in parts) {
      if (is.numeric(data[[p]])) numeric_vars <- union(numeric_vars, p)
      if (is.factor(data[[p]]))  factor_vars  <- union(factor_vars, p)
    }
  }
  
  list(numeric_vars = numeric_vars, factor_vars = factor_vars)
}

# ------------------------------------------------------------
# Helper: Build prediction grid
# ------------------------------------------------------------
build_prediction_grid <- function(data, numeric_vars, factor_vars, numeric_setting) {
  
  # numeric values
  nv <- lapply(numeric_vars, function(v) {
    if (numeric_setting == "mean") mean(data[[v]], na.rm = TRUE)
    else median(data[[v]], na.rm = TRUE)
  })
  names(nv) <- numeric_vars
  
  # factor level combinations
  if (length(factor_vars) > 0) {
    fac_levels <- lapply(factor_vars, function(v) levels(data[[v]]))
    names(fac_levels) <- factor_vars
    base_grid <- expand.grid(fac_levels, stringsAsFactors = FALSE)
  } else {
    base_grid <- data.frame(dummy = 1)[FALSE, ]
    base_grid <- data.frame(row = 1)
    base_grid$row <- NULL
  }
  
  # add numeric values
  for (v in numeric_vars) base_grid[[v]] <- nv[[v]]
  
  base_grid$numeric_setting <- numeric_setting
  
  base_grid
}

# ============================================================
# MAIN LOOP OVER DATASETS
# ============================================================
for (i in seq_along(datasets)) {
  
  ds_name <- names(datasets)[i]
  df <- datasets[[i]]
  
  cat("\n=== Predictions for dataset:", ds_name, "===\n")
  
  uniq_map <- model_uniqueness[[ds_name]]
  if (is.null(uniq_map)) next
  
  pred_all <- data.frame()
  
  # Loop over selection methods
  for (method in method_order) {
    
    # skip duplicate models
    if (uniq_map[[method]] != method) next
    
    method_res <- selected_models[[ds_name]][[method]]
    if (is.null(method_res)) next
    
    form <- get_model_formula(method_res)
    if (is.null(form)) next
    
    fit <- tryCatch(lm(form, data = df), error = \(e) NULL)
    if (is.null(fit)) next
    
    comps <- get_model_components(fit, df)
    
    numeric_vars <- comps$numeric_vars
    factor_vars  <- comps$factor_vars
    
    # Build mean & median prediction grids
    grid <- bind_rows(
      build_prediction_grid(df, numeric_vars, factor_vars, "mean"),
      build_prediction_grid(df, numeric_vars, factor_vars, "median")
    )
    grid$Method <- method
    
    # restrict newdata to model variables
    mf <- model.frame(fit)
    vars_in_model <- names(mf)[names(mf) != "risk"]
    newdata <- grid[, intersect(vars_in_model, names(grid)), drop = FALSE]
    
    # compute intervals
    pi <- predict(fit, newdata = newdata, interval = "prediction")
    ci <- predict(fit, newdata = newdata, interval = "confidence")
    
    # Combine into interval notation
    interval_PI <- paste0("[", round(pi[, "lwr"], 3), ", ", round(pi[, "upr"], 3), "]")
    interval_CI <- paste0("[", round(ci[, "lwr"], 3), ", ", round(ci[, "upr"], 3), "]")
    
    # Build output row for this model
    out <- grid %>%
      select(Method, numeric_setting, all_of(factor_vars)) %>%
      mutate(
        Prediction_Interval = interval_PI,
        Confidence_Interval = interval_CI
      )
    
    pred_all <- bind_rows(pred_all, out)
  }
  
  # ----------------------------------------------------------
  # WRITE OUTPUT FILE
  # ----------------------------------------------------------
  ds_dir <- file.path(pred_path, ds_name)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  
  outfile <- file.path(ds_dir, paste0(ds_name, "_predictions.csv"))
  write.csv(pred_all, outfile, row.names = FALSE)
  
  cat("✔ Saved:", outfile, "\n")
}
