# Collinearity.R
# Step: Collinearity diagnostics and eigprop-based variable reduction

library(tidyverse)
library(car)
library(mctest)
library(GGally)

# Global ggplot theme (optional)
theme_set(
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85"),
      plot.margin = margin(10, 10, 10, 10)
    )
)

# --- Paths ---
base_path <- "/cloud/project/Fall 2025/Final Project Fall 2025"
data_path <- file.path(base_path, "Data")
results_path <- file.path(base_path, "Results", "Collinearity")
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------
# Helper: load a single dataset
# ---------------------------------------------
load_dataset <- function(f) {
  e <- new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]
}

# ---------------------------------------------
# Helper: drop single-level factor variables
# ---------------------------------------------
drop_single_level_factors <- function(df) {
  fvars <- names(Filter(is.factor, df))
  single <- fvars[sapply(df[fvars], function(x) length(unique(x)) < 2)]
  if (length(single) > 0) {
    message("Dropping single-level factors: ", paste(single, collapse = ", "))
    df <- df |> select(-all_of(single))
  }
  df
}

# ---------------------------------------------
# Helper: remove aliased TERMS (not just columns)
# Returns the names of terms to drop, e.g. "gender", "age:gender"
# ---------------------------------------------
aliased_terms_to_drop <- function(formula_full, df) {
  # Fit full model
  fit <- lm(formula_full, data = df)
  
  # Design matrix and mapping column -> term index
  mm     <- model.matrix(fit)
  assign <- attr(mm, "assign")              # integer term indices (0 = intercept)
  terms  <- attr(terms(fit), "term.labels") # names of model terms
  
  # Alias info: which columns are linearly dependent
  comp <- alias(fit)$Complete
  if (is.null(comp)) return(character())
  
  aliased_cols <- rownames(comp)  # these are column names of mm
  if (is.null(aliased_cols) || length(aliased_cols) == 0) return(character())
  
  # Map each aliased column back to its term
  term_ids <- unique(assign[colnames(mm) %in% aliased_cols])
  term_ids <- term_ids[term_ids > 0]  # drop intercept index 0, just in case
  
  if (length(term_ids) == 0) return(character())
  
  unique(terms[term_ids])
}

# --- Thresholds ---
ci_thresh <- 30
pi_thresh <- 0.50

# Initialize output file ONCE before the loop (outside for-loop)
remaining_summary_file <- file.path(results_path, "Remaining_Variables_All_Datasets.txt")
if (file.exists(remaining_summary_file)) file.remove(remaining_summary_file)

# ---------------------------------------------
# MAIN LOOP
# ---------------------------------------------
files <- list.files(data_path, pattern = "\\.RData$", full.names = TRUE)

for (file in files) {
  
  dataset_name <- tools::file_path_sans_ext(basename(file))
  message("\n===================================================")
  message("Processing dataset: ", dataset_name)
  message("===================================================\n")
  
  # Output folder
  out_dir <- file.path(results_path, dataset_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------
  # Load dataset
  # -------------------------------------------
  df <- load_dataset(file)
  
  if (!"risk" %in% names(df)) {
    warning("Dataset missing 'risk': ", dataset_name)
    next
  }
  
  # -------------------------------------------
  # Remove single-level factors
  # -------------------------------------------
  df <- drop_single_level_factors(df)
  
  # -------------------------------------------
  # Identify predictor variables
  # -------------------------------------------
  factor_vars  <- names(Filter(is.factor, df))
  numeric_vars <- setdiff(names(Filter(is.numeric, df)), "risk")
  
  # -------------------------------------------
  # Full predictor list
  # -------------------------------------------
  all_vars <- c(numeric_vars, factor_vars)
  
  # -------------------------------------------
  # Build full formula
  # -------------------------------------------
  full_formula <- as.formula(
    paste("risk ~", paste(all_vars, collapse = " + "))
  )
  
  # -------------------------------------------
  # Remove aliased TERMS BEFORE eigprop
  # -------------------------------------------
  terms_to_drop <- aliased_terms_to_drop(full_formula, df)
  
  if (length(terms_to_drop) > 0) {
    message("Dropping aliased terms: ", paste(terms_to_drop, collapse = ", "))
    
    # all_vars contains term labels: numeric vars, factor vars, and interactions
    all_vars     <- setdiff(all_vars, terms_to_drop)
    numeric_vars <- setdiff(numeric_vars, terms_to_drop)
    factor_vars  <- setdiff(factor_vars, terms_to_drop)
    
    if (length(all_vars) == 0) {
      warning("All predictors dropped due to aliasing in ", dataset_name)
    } else {
      full_formula <- as.formula(
        paste("risk ~", paste(all_vars, collapse = " + "))
      )
    }
  }
  
  # -------------------------------------------
  # EIGPROP LOOP (INTERACTIVE)
  # -------------------------------------------
  remaining <- all_vars
  removed   <- character()
  iter      <- 1
  
  eig_log <- data.frame(
    Iteration = integer(),
    Variable_Removed = character(),
    CI_Value = numeric(),
    PI_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  repeat {
    
    message("\n--------------------------------------------")
    message("Iteration ", iter, " for dataset ", dataset_name)
    message("--------------------------------------------")
    
    # Refit model
    f_curr <- as.formula(
      paste("risk ~", paste(remaining, collapse = " + "))
    )
    lmod_curr <- lm(f_curr, df)
    
    # Compute eigprop
    eig <- eigprop(lmod_curr)
    
    # =====================================================
    # SAFELY DEFINE ROW_INDEX AFTER EACH EIGPROP CALL
    # =====================================================
    if (is.null(eig) || is.null(eig$ci)) {
      message("eigprop returned NULL — stopping collinearity process for this dataset.")
      break
    }
    
    # Number of rows in eigprop output
    nrows <- length(eig$ci)
    
    if (nrows == 0) {
      message("eigprop produced zero rows — stopping collinearity process for this dataset.")
      break
    }
    
    # Start on the LAST ROW (highest CI)
    row_index <- nrows
    
    # Start with highest CI row
    order_rows <- order(ci_vals, decreasing = TRUE)
    dropped_this_iter <- FALSE
    
    # =============================================================
    # Convert column-level PI values from eigprop → TERM-level PI
    # =============================================================
    terms_all  <- attr(lmod_curr$terms, "term.labels") # character term labels
    
    # map variable -> model columns
    var_to_cols <- lapply(terms_all, function(v) {
      cols <- grep(paste0("^", v), colnames(eig$pi), value = TRUE)
      cols
    })
    names(var_to_cols) <- terms_all
    
    var_tstat  <- c()
    var_pvalue <- c()
    
    # Get full model summary once
    sm <- summary(lmod_curr)$coefficients   # matrix: term × (Estimate, StdErr, t, p)
    
    # Prevent undefined or invalid row_index from causing errors
    if (is.null(row_index) || !is.numeric(row_index)) break
    if (row_index > length(eig$ci)) row_index <- length(eig$ci)
    
    while (row_index >= 1 && eig$ci[row_index] >= ci_thresh) {
      # Compute group PI
      group_PI <- sapply(names(var_to_cols), function(v) {
        cols <- var_to_cols[[v]]
        sum(eig$pi[row_index, cols], na.rm = TRUE)
      })
      
      # Identify all terms exceeding the PI threshold
      terms_above_thresh <- names(group_PI[group_PI > pi_thresh])
      
      # =============================================================
      # Choose TERM to drop based on TERM-level PI values
      # =============================================================
      
      if (length(terms_above_thresh) == 0) {
        # No terms exceed the PI threshold → move up a row
        row_index <- row_index - 1
        next
      }
      
      # =============================================================
      # Display TERMS above PI threshold *with PI values*
      # =============================================================
      if (length(terms_above_thresh) > 1) {
        
        for (var in terms_above_thresh) {
          
          # all columns in design matrix associated with this variable
          cols <- var_to_cols[[var]]
          
          # Find coefficient rows in model summary corresponding to those columns
          # Matching by prefix handles dummy variables like incomeMedium, incomeHigh
          coef_rows <- do.call(c, lapply(cols, function(colname) {
            grep(paste0("^", colname, "$"), rownames(sm))
          }))
          
          # Remove NAs and proceed
          coef_rows <- coef_rows[!is.na(coef_rows)]
          
          if (length(coef_rows) == 0) {
            # If no coefficient rows found (rare), set NA
            var_tstat[var]  <- NA
            var_pvalue[var] <- NA
            next
          }
          
          # Extract t‐statistics and p‐values for all dummy columns
          t_stats <- abs(sm[coef_rows, "t value"])
          p_vals  <- sm[coef_rows, "Pr(>|t|)"]
          
          # Aggregate per variable:
          #   - t-statistic = max absolute t across levels
          #   - p-value     = min across levels
          var_tstat[var]  <- max(t_stats, na.rm = TRUE)
          var_pvalue[var] <- min(p_vals, na.rm = TRUE)
        }
        
        # =============================================================
        # DISPLAY TABLE (PI + t-stat + p-value)
        # =============================================================
        
        table_to_print <- data.frame(
          Variable   = terms_above_thresh,
          Group_PI   = round(group_PI[terms_above_thresh], 4),
          T_stat     = round(var_tstat[terms_above_thresh], 4),
          P_value    = signif(var_pvalue[terms_above_thresh], 4)
        )
        
        message("Multiple TERMS exceed PI threshold (", pi_thresh, "):")
        
        print(table_to_print, row.names = FALSE)
        
        drop_term <- readline("Enter TERM to drop: ")
        
        if (!(drop_term %in% terms_above_thresh)) {
          message("Invalid choice; dropping worst TERM by PI value: ", names(pi_sub)[1])
          drop_term <- names(pi_sub)[1]
        }
        
      } else {
        
        # Only one TERM exceeds threshold
        drop_term <- terms_above_thresh[1]
        message("Dropping TERM: ", drop_term, 
                " (PI = ", round(pi_by_term[drop_term], 4), ")")
      }
      
      eig_log <- rbind(
        eig_log,
        data.frame(
          Iteration = iter,
          Variable_Removed = drop_term,
          CI_Value = ci_vals[r],
          PI_Value = pi_row[drop_term]
        )
      )
      
      message("Dropping TERM: ", drop_term)
      
      # Remove entire term (works for numeric, factor, interaction)
      remaining <- setdiff(remaining, drop_term)
      
      # Restart eigprop loop after dropping term
      dropped_this_iter <- TRUE
      break
    }
    
    if (!dropped_this_iter) {
      message("No more variables meet CI/PI thresholds.")
      break
    }
    
    if (length(remaining) == 0) break
    
    iter <- iter + 1
  }
  
  # -------------------------------------------
  # Save eigprop log
  # -------------------------------------------
  write.csv(
    eig_log,
    file.path(out_dir, paste0("eigprop_log_", dataset_name, ".csv")),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------------------------
  # Append remaining variables for this dataset to the master summary file
  # ------------------------------------------------------------------------------
  
  cat("\n============================================================\n",
      "Dataset: ", dataset_name, "\n",
      "============================================================\n",
      file = remaining_summary_file, append = TRUE)
  
  if (length(remaining) == 0) {
    cat("NO REMAINING VARIABLES (all dropped due to collinearity)\n",
        file = remaining_summary_file, append = TRUE)
  } else {
    cat("Remaining Variables (sorted):\n", file = remaining_summary_file, append = TRUE)
    cat(paste0(" - ", sort(remaining), collapse = "\n"),
        "\n\n",
        file = remaining_summary_file,
        append = TRUE)
  }
  
  # ============================================================
  #           POST-EIGPROP ANALYSIS  (add this block)
  # ============================================================
  
  if (length(remaining) > 0) {
    
    
    final_lm <- lm(final_formula, data = df)
    
    # -------------------------------
    # VIF values for final model
    # -------------------------------
    vif_out <- tryCatch(vif(final_lm), error = function(e) NULL)
    if (!is.null(vif_out)) {
      vif_df <- if (is.matrix(vif_out)) {
        data.frame(Variable = rownames(vif_out), VIF = vif_out[,1])
      } else {
        data.frame(Variable = names(vif_out), VIF = as.numeric(vif_out))
      }
      write.csv(
        vif_df,
        file.path(out_dir, paste0("vif_", dataset_name, ".csv")),
        row.names = FALSE
      )
    }
    
    # -------------------------------
    # Remaining numeric variables
    # -------------------------------
    remaining_numeric <- intersect(names(Filter(is.numeric, df)), remaining)
    
    # -------------------------------
    # Correlation matrix
    # -------------------------------
    if (length(remaining_numeric) > 1) {
      cor_mat <- cor(df[, remaining_numeric, drop = FALSE], use = "pairwise.complete.obs")
      write.csv(
        cor_mat,
        file.path(out_dir, paste0("correlation_matrix_", dataset_name, ".csv"))
      )
    }
    
    # -------------------------------
    # Pairs plot (GGally::ggpairs)
    # -------------------------------
    if (length(remaining_numeric) > 1) {
      p <- GGally::ggpairs(df[, remaining_numeric, drop = FALSE])
      ggsave(
        file.path(out_dir, paste0("pairs_plot_", dataset_name, ".png")),
        p,
        width = 10, height = 10, dpi = 300
      )
    }
    
    # -------------------------------
    # Factor–factor heatmaps
    # -------------------------------
    remaining_factor <- intersect(names(Filter(is.factor, df)), remaining)
    
    if (length(remaining_factor) > 1) {
      for (i in 1:(length(remaining_factor)-1)) {
        for (j in (i+1):length(remaining_factor)) {
          
          f1 <- remaining_factor[i]
          f2 <- remaining_factor[j]
          
          tab <- table(df[[f1]], df[[f2]])
          df_tab <- as.data.frame(tab)
          names(df_tab) <- c("Var1", "Var2", "Freq")
          
          g <- ggplot(df_tab, aes(x = Var1, y = Var2, fill = Freq)) +
            geom_tile() +
            geom_text(aes(label = Freq), color = "white") +
            scale_fill_gradient() +
            labs(
              title = paste("Factor–Factor Heatmap:", f1, "vs", f2),
              x = f1,
              y = f2
            ) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          
          ggsave(
            file.path(out_dir, paste0("heatmap_", f1, "_vs_", f2, ".png")),
            g,
            width = 6, height = 5, dpi = 300
          )
        }
      }
    }
  }
  
  
  
  cat("Collinearity analysis completed for dataset:", dataset_name, "\n")
}
