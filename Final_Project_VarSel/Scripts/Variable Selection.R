# ==========================================
# Variable Selection - After Collinearity
# ==========================================

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("leaps")) install.packages("leaps")
if (!require("ggplot2")) install.packages("ggplot2")

library(tidyverse)
library(leaps)
library(ggplot2)

# ======================================================
# Read dataset name from command-line argument (SLURM)
# ======================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No dataset name supplied. Usage: Rscript Variable_Selection.R <DatasetName>")
}

dataset_to_run <- args[1]
message("Running variable selection for dataset: ", dataset_to_run)

# ---- Paths ----
base_path    <- "~/Math 5387 Regression Analysis/Final Project Fall 2025"
data_path    <- file.path(base_path, "Data")
collin_path  <- file.path(base_path, "Results", "Collinearity")
vs_path      <- file.path(base_path, "Results", "Variable_Selection")
dir.create(vs_path, recursive = TRUE, showWarnings = FALSE)

# ======================================================
# Load ONLY the dataset passed from SLURM
# ======================================================
dataset_file <- file.path(data_path, paste0(dataset_to_run, ".RData"))

if (!file.exists(dataset_file)) {
  stop("Dataset file not found: ", dataset_file)
}

env <- new.env()
load(dataset_file, envir = env)
df <- env[[ls(env)[1]]]  # extract dataset

variable_selection <- function(data,
                               remaining_vars,
                               response_var = "risk") {
  # keep only vars that actually exist
  remaining_vars <- remaining_vars[remaining_vars %in% names(data)]
  
  # --------------------------
  # Split variable types
  # --------------------------
  interaction_vars <- remaining_vars[grepl(":", remaining_vars)]
  main_vars        <- setdiff(remaining_vars, interaction_vars)
  
  numeric_main <- main_vars[sapply(data[main_vars], is.numeric)]
  factor_main  <- main_vars[sapply(data[main_vars], is.factor)]
  
  # all terms in full model
  all_terms <- c(numeric_main, factor_main, interaction_vars)
  
  if (length(all_terms) == 0) {
    stop("No remaining predictors provided to variable_selection().")
  }
  
  # full formula: risk ~ all main + interaction terms
  formula_full <- as.formula(
    paste(response_var, "~", paste(all_terms, collapse = " + "))
  )
  
  # null model (intercept only)
  formula_null <- as.formula(paste(response_var, "~ 1"))
  
  # fit full & null
  fit_full <- lm(formula_full, data = data)
  fit_null <- lm(formula_null, data = data)
  
  n <- nrow(data)
  
  results <- list()
  
  # Helper: extract variable + order + criterion
  get_p_table <- function(fit) {
    s <- summary(fit)
    coefs <- s$coefficients
    if (nrow(coefs) <= 1) {
      return(data.frame(variable = character(),
                        order = integer(),
                        criterion_value = numeric()))
    }
    pvals <- coefs[-1, 4]
    vars  <- rownames(coefs)[-1]
    data.frame(
      variable = vars,
      order = seq_along(vars),
      criterion_value = as.numeric(pvals),
      stringsAsFactors = FALSE
    )
  }
  
  get_AIC_table <- function(fit) {
    a <- AIC(fit)
    vars <- names(coef(fit))[-1]
    data.frame(
      variable = vars,
      order = seq_along(vars),
      criterion_value = rep(a, length(vars)),
      stringsAsFactors = FALSE
    )
  }
  
  # --------------------------
  # Stepwise-like methods via step()
  # (Note: step uses AIC, we preserve your naming)
  # --------------------------
  
  # backward "p" (your original label; actually BIC-ish with k=log(n))
  backward_p <- step(fit_full, direction = "backward", k = log(n), trace = 0)
  results$backward_p <- list(
    df = get_p_table(backward_p),
    final_model = formula(backward_p)
  )
  
  # forward "p"
  forward_p <- step(fit_null, scope = formula_full,
                    direction = "forward", k = log(n), trace = 0)
  results$forward_p <- list(
    df = get_p_table(forward_p),
    final_model = formula(forward_p)
  )
  
  # stepwise "p"
  stepwise_p <- step(fit_null, scope = formula_full,
                     direction = "both", k = log(n), trace = 0)
  results$stepwise_p <- list(
    df = get_p_table(stepwise_p),
    final_model = formula(stepwise_p)
  )
  
  # backward AIC
  backward_AIC <- step(fit_full, direction = "backward", trace = 0)
  results$backward_AIC <- list(
    df = get_AIC_table(backward_AIC),
    final_model = formula(backward_AIC)
  )
  
  # forward AIC
  forward_AIC <- step(fit_null, scope = formula_full,
                      direction = "forward", trace = 0)
  results$forward_AIC <- list(
    df = get_AIC_table(forward_AIC),
    final_model = formula(forward_AIC)
  )
  
  # stepwise AIC
  stepwise_AIC <- step(fit_null, scope = formula_full,
                       direction = "both", trace = 0)
  results$stepwise_AIC <- list(
    df = get_AIC_table(stepwise_AIC),
    final_model = formula(stepwise_AIC)
  )
  
  # --------------------------
  # Subset selection via regsubsets
  # with dummy-expanded design matrix
  # --------------------------
  
  # Design matrix for full model (no intercept)
  mm <- model.matrix(formula_full, data = data)
  # drop intercept column if present
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  
  y <- data[[response_var]]
  
  reg_fit <- regsubsets(
    x = mm,
    y = y,
    nbest = 1,
    nvmax = ncol(mm)
  )
  reg_summary <- summary(reg_fit)
  
  # Mapping columns → original term labels
  assign_vec   <- attr(mm, "assign")     # term index per column
  terms_full   <- terms(formula_full)
  term_labels  <- attr(terms_full, "term.labels")
  
  # Number of predictors (including intercept) per model row
  p_vec <- apply(reg_summary$which, 1, sum)
  rss   <- reg_summary$rss
  nmod  <- length(rss)
  
  aic_vec  <- n * log(rss / n) + 2 * p_vec
  bic_vec  <- reg_summary$bic
  adjr_vec <- reg_summary$adjr2
  rmse_vec <- sqrt(rss / n)
  cp_vec   <- reg_summary$cp
  
  # helper: given index for best model, return unique term labels
  terms_from_model_idx <- function(idx) {
    which_row <- reg_summary$which[idx, ]
    # drop intercept column (assumed first)
    if ("(Intercept)" %in% colnames(reg_summary$which)) {
      which_row <- which_row[colnames(reg_summary$which) != "(Intercept)"]
    }
    col_idx <- which(which_row)
    if (length(col_idx) == 0) return(character())
    # map each column to term label via assign
    unique(term_labels[assign_vec[col_idx]])
  }
  
  # Best models indices
  idx_AIC  <- which.min(aic_vec)
  idx_BIC  <- which.min(bic_vec)
  idx_AdjR2 <- which.max(adjr_vec)
  idx_RMSE <- which.min(rmse_vec)
  idx_Cp   <- which.min(cp_vec)
  
  # Build data.frames for subsets methods
  make_subset_df <- function(idx, crit_vec, label) {
    vars <- terms_from_model_idx(idx)
    if (length(vars) == 0) {
      return(data.frame(variable = character(),
                        order = integer(),
                        criterion_value = numeric()))
    }
    data.frame(
      variable        = vars,
      order           = seq_along(vars),
      criterion_value = rep(crit_vec[idx], length(vars)),
      stringsAsFactors = FALSE
    )
  }
  
  results$subsets_AIC <- list(
    df = make_subset_df(idx_AIC, aic_vec, "AIC"),
    final_model = NULL
  )
  results$subsets_BIC <- list(
    df = make_subset_df(idx_BIC, bic_vec, "BIC"),
    final_model = NULL
  )
  results$subsets_AdjR2 <- list(
    df = make_subset_df(idx_AdjR2, adjr_vec, "AdjR2"),
    final_model = NULL
  )
  results$subsets_RMSE <- list(
    df = make_subset_df(idx_RMSE, rmse_vec, "RMSE"),
    final_model = NULL
  )
  results$subsets_Cp <- list(
    df = make_subset_df(idx_Cp, cp_vec, "Cp"),
    final_model = NULL
  )
  
  # Also attach reg_summary & criteria for plotting
  results$reg_summary <- list(
    which = reg_summary$which,
    AIC   = aic_vec,
    BIC   = bic_vec,
    AdjR2 = adjr_vec,
    RMSE  = rmse_vec,
    Cp    = cp_vec
  )
  
  return(results)
}

plot_subset_diagnostics <- function(reg_summary_list, out_dir, dataset_name) {
  which_mat <- reg_summary_list$which
  # number of predictors (incl. intercept)
  p_vec <- apply(which_mat, 1, sum)
  
  df_plot <- data.frame(
    ModelSize = p_vec,
    AIC  = reg_summary_list$AIC,
    BIC  = reg_summary_list$BIC,
    AdjR2 = reg_summary_list$AdjR2,
    RMSE = reg_summary_list$RMSE,
    Cp   = reg_summary_list$Cp
  )
  
  # AIC
  g_aic <- ggplot(df_plot, aes(x = ModelSize, y = AIC)) +
    geom_line() + geom_point() +
    labs(title = paste("AIC vs Model Size –", dataset_name),
         x = "Number of Predictors", y = "AIC") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "subsets_AIC_vs_size.png"),
         g_aic, width = 6, height = 4, dpi = 300)
  
  # BIC
  g_bic <- ggplot(df_plot, aes(x = ModelSize, y = BIC)) +
    geom_line() + geom_point() +
    labs(title = paste("BIC vs Model Size –", dataset_name),
         x = "Number of Predictors", y = "BIC") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "subsets_BIC_vs_size.png"),
         g_bic, width = 6, height = 4, dpi = 300)
  
  # AdjR2
  g_adj <- ggplot(df_plot, aes(x = ModelSize, y = AdjR2)) +
    geom_line() + geom_point() +
    labs(title = paste("Adjusted R² vs Model Size –", dataset_name),
         x = "Number of Predictors", y = "Adjusted R²") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "subsets_AdjR2_vs_size.png"),
         g_adj, width = 6, height = 4, dpi = 300)
  
  # RMSE
  g_rmse <- ggplot(df_plot, aes(x = ModelSize, y = RMSE)) +
    geom_line() + geom_point() +
    labs(title = paste("RMSE vs Model Size –", dataset_name),
         x = "Number of Predictors", y = "RMSE") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "subsets_RMSE_vs_size.png"),
         g_rmse, width = 6, height = 4, dpi = 300)
  
  # Cp
  g_cp <- ggplot(df_plot, aes(x = ModelSize, y = Cp)) +
    geom_line() + geom_point() +
    labs(title = paste("Mallows Cp vs Model Size –", dataset_name),
         x = "Number of Predictors", y = "Cp") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "subsets_Cp_vs_size.png"),
         g_cp, width = 6, height = 4, dpi = 300)
}

write_selected_variables_per_dataset <- function(
    selected_vars_list,
    dataset_names,
    output_path,
    combined_model_file = "Final_Models_All_Datasets.txt"
) {
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  final_model_order <- c(
    "backward_p", "forward_p", "stepwise_p",
    "backward_AIC", "forward_AIC", "stepwise_AIC",
    "subsets_AIC", "subsets_BIC", "subsets_AdjR2",
    "subsets_RMSE", "subsets_Cp"
  )
  
  pretty_names <- list(
    backward_p   = "Backward_pvalue",
    forward_p    = "Forward_pvalue",
    stepwise_p   = "Stepwise_pvalue",
    backward_AIC = "Backward_AIC",
    forward_AIC  = "Forward_AIC",
    stepwise_AIC = "Stepwise_AIC",
    subsets_AIC  = "Subsets_AIC",
    subsets_BIC  = "Subsets_BIC",
    subsets_AdjR2 = "Subsets_AdjR2",
    subsets_RMSE = "Subsets_RMSE",
    subsets_Cp   = "Subsets_Cp"
  )
  
  combined_file_path <- file.path(output_path, combined_model_file)
  sink(combined_file_path)
  cat("Final Selected Models for All Datasets\n\n")
  
  for (i in seq_along(selected_vars_list)) {
    ds_name <- dataset_names[i]
    methods <- selected_vars_list[[i]]
    
    # Collect all variables across methods
    all_vars <- unique(unlist(lapply(methods[final_model_order], function(m) {
      if (is.null(m) || is.null(m$df)) return(NULL)
      m$df$variable
    })))
    out <- data.frame(variable = all_vars, stringsAsFactors = FALSE)
    
    for (method_key in final_model_order) {
      pretty <- pretty_names[[method_key]]
      m <- methods[[method_key]]
      
      if (!is.null(m) && !is.null(m$df) &&
          all(c("variable","order","criterion_value") %in% names(m$df))) {
        
        dfm <- m$df[order(m$df$order), ]
        
        tmp_order <- dfm[, c("variable", "order")]
        names(tmp_order)[2] <- paste0(pretty, "_Order")
        out <- merge(out, tmp_order, by = "variable", all.x = TRUE)
        
        tmp_val <- dfm[, c("variable", "criterion_value")]
        names(tmp_val)[2] <- paste0(pretty, "_Value")
        out <- merge(out, tmp_val, by = "variable", all.x = TRUE)
      } else {
        out[, paste0(pretty, "_Order")] <- NA
        out[, paste0(pretty, "_Value")] <- NA
      }
    }
    
    # Write CSV with two header rows
    outfile <- file.path(output_path, paste0(ds_name, "_Selected_Variables.csv"))
    col_names <- names(out)
    first_header  <- col_names
    second_header <- col_names
    second_header[1] <- ""
    for (j in 2:length(col_names)) {
      if (grepl("_Order$", col_names[j])) second_header[j] <- "Order"
      if (grepl("_Value$", col_names[j])) second_header[j] <- "Value"
    }
    con <- file(outfile, "w")
    writeLines(paste(first_header, collapse = ","), con)
    writeLines(paste(second_header, collapse = ","), con)
    write.table(out, con, sep = ",", row.names = FALSE,
                col.names = FALSE, na = "")
    close(con)
    message("Saved variable table: ", outfile)
    
    # Append final model formulas
    cat("===== Dataset: ", ds_name, " =====\n\n", sep = "")
    for (method_name in final_model_order) {
      m <- methods[[method_name]]
      if (!is.null(m) && !is.null(m$final_model)) {
        pretty <- pretty_names[[method_name]]
        cat("---- ", pretty, " ----\n", sep = "")
        cat(deparse(m$final_model), "\n\n")
      }
    }
  }
  
  sink()
  message("Saved combined final models: ", combined_file_path)
}

identify_unique_models <- function(selected_models_for_dataset, method_order) {
  
  normalize_term <- function(term) {
    # If interaction term (contains ":")
    if (grepl(":", term)) {
      parts <- unlist(strsplit(term, ":"))
      parts <- sort(parts)   # standardize A:B to always alphabetical
      return(paste(parts, collapse = ":"))
    } else {
      return(term)
    }
  }
  
  extract_vars <- function(model_result) {
    # Case 1: result has full formula
    if (!is.null(model_result$final_model)) {
      vars <- attr(terms(model_result$final_model), "term.labels")
    }
    # Case 2: subsets results (stored in a df)
    else if (!is.null(model_result$df)) {
      vars <- model_result$df$variable
    } else {
      return(NULL)
    }
    
    # Normalize interaction terms, remove whitespace
    vars <- gsub(" ", "", vars)
    vars <- sapply(vars, normalize_term)
    
    # Sort → order independent
    vars <- sort(vars)
    vars
  }
  
  model_keys   <- list()
  canonical    <- list()
  unique_models <- list()
  
  for (method in method_order) {
    result <- selected_models_for_dataset[[method]]
    if (is.null(result)) next
    
    vars <- extract_vars(result)
    if (is.null(vars)) next
    
    # Build a stable, order-independent key
    key <- paste(vars, collapse = " + ")
    
    model_keys[[method]] <- key
  }
  
  seen <- character(0)
  
  for (method in method_order) {
    key <- model_keys[[method]]
    if (is.null(key)) next
    
    if (!(key %in% seen)) {
      # First time this model appears — record representative method name
      unique_models[[key]] <- method
      seen <- c(seen, key)
    }
    
    # Map method → canonical name
    canonical[[method]] <- unique_models[[key]]
  }
  
  list(
    model_keys     = model_keys,     # raw keys for each method
    canonical      = canonical,      # mapping method → representative method
    unique_models  = unique_models   # mapping key → representative method
  )
}

selected_vars_list <- vector("list", length(datasets))
names(selected_vars_list) <- names(datasets)


# Collinearity result directory for this dataset
ds_collin_dir <- file.path(collin_path, dataset_to_run)
removed_file  <- file.path(ds_collin_dir, "eigprop_removed_variables.csv")
  
# Infer remaining vars: all non-response predictors minus removed
all_predictors <- setdiff(names(df), "risk")
  
removed <- if (file.exists(removed_file)) {
  read.csv(removed_file, stringsAsFactors = FALSE)$Variable_Removed
} else {
  character(0)
}

remaining_vars <- setdiff(all_predictors, removed)
  
sel <- variable_selection(
  data = df,
  remaining_vars = remaining_vars,
  response_var = "risk"
)

# make per-dataset directory for VS
ds_vs_dir <- file.path(vs_path, dataset_to_run)
dir.create(ds_vs_dir, recursive = TRUE, showWarnings = FALSE)
  
# subset diagnostic plots
plot_subset_diagnostics(sel$reg_summary, ds_vs_dir, ds_name)

selected_vars_list <- list(sel)
names(selected_vars_list) <- dataset_to_run

# Write tables + combined final model formulas
write_selected_variables_per_dataset(
  selected_vars_list = selected_vars_list,
  dataset_names = dataset_to_run,
  output_path = vs_path,
  combined_model_file = "Final_Models_All_Datasets.txt"
)

# Create lookup table
model_uniqueness <- list()

for (i in seq_along(selected_vars_list)) {
  ds_name <- names(selected_vars_list)[i]
  model_uniqueness[[ds_name]] <- identify_unique_models(
    selected_models_for_dataset = selected_vars_list[[i]],
    method_order = c(
      "backward_p", "forward_p", "stepwise_p",
      "backward_AIC", "forward_AIC", "stepwise_AIC",
      "subsets_AIC", "subsets_BIC", "subsets_AdjR2",
      "subsets_RMSE", "subsets_Cp"
    )
  )
}

saveRDS(model_uniqueness,
        file = file.path(results_path, "model_uniqueness_map.rds"))
