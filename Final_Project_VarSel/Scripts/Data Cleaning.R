# --------------------------
# Data Cleaning - Regression Analysis
# --------------------------

## Packages
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# --- Paths ---
base_path <- "/cloud/project/Fall 2025/Final Project Fall 2025"
data_path <- file.path(base_path,"Data")
data_file <- file.path(data_path, "diabetes_dataset.csv")
results_path <- file.path(base_path, "Results")
summaries_path <- file.path(results_path, "Summaries")
dir.create(summaries_path, showWarnings = FALSE)

# --- Load Raw Data ---
diabetes <- read.csv(data_file, header = TRUE)

# --- Rename & clean variables ---
diabetes_clean <- diabetes |>
  rename(
    education = education_level,
    income = income_level,
    employment = employment_status,
    smoking = smoking_status,
    alcohol_use = alcohol_consumption_per_week,
    physical_activity = physical_activity_minutes_per_week,
    sleep = sleep_hours_per_day,
    screen_time = screen_time_hours_per_day,
    diabetes_hist = family_history_diabetes,
    hypertension_hist = hypertension_history,
    cardio_hist = cardiovascular_history,
    waist_hip = waist_to_hip_ratio,
    hr = heart_rate,
    total_chol = cholesterol_total,
    hdl = hdl_cholesterol,
    ldl = ldl_cholesterol,
    insulin = insulin_level,
    risk = diabetes_risk_score,
    type = diabetes_stage,
    diagnosis = diagnosed_diabetes
  )

# --- Convert factor variables ---
factor_vars <- c("gender","ethnicity","education","income","employment",
                 "smoking","diabetes_hist","hypertension_hist",
                 "cardio_hist","type","diagnosis")

diabetes_clean <- diabetes_clean |>
  mutate(across(all_of(factor_vars), as.factor))

# --- Rename levels for binary factors ---
levels(diabetes_clean$diabetes_hist) = c("No","Yes")
levels(diabetes_clean$hypertension_hist) = c("No","Yes")
levels(diabetes_clean$cardio_hist) = c("No","Yes")
levels(diabetes_clean$diagnosis) = c("Undiagnosed","Diagnosed")

# ---------------------------
# Split by type and diagnosis
# ---------------------------
type_levels <- levels(diabetes_clean$type)
diagnosis_levels <- levels(diabetes_clean$diagnosis)
combinations <- expand.grid(type = type_levels, diagnosis = diagnosis_levels)

# Initialize counts dataframe
counts_df <- data.frame(
  type_diagnosis = character(),
  count = integer(),
  stringsAsFactors = FALSE
)

# Function to summarize numeric variables
summarize_numeric <- function(x) {
  c(
    Min = min(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Mean = mean(x, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE),
    Max = max(x, na.rm = TRUE)
  )
}

# Keep track of skipped datasets
skipped_combinations <- c()

for(i in 1:nrow(combinations)) {
  current_type <- combinations$type[i]
  current_diag <- combinations$diagnosis[i]
  
  current_data <- diabetes_clean |>
    filter(type == current_type & diagnosis == current_diag)
  
  # Skip if no observations
  if(nrow(current_data) == 0) {
    skipped_combinations <- c(skipped_combinations, paste(current_type, current_diag, sep = "_"))
    next
  }
  
  # Drop type and diagnosis (already done, but ensure)
  current_data <- current_data |>
    select(-type, -diagnosis, -total_chol)
  
  # Drop factor variables with only one level
  factor_cols <- names(select(current_data, where(is.factor)))
  single_value_factors <- factor_cols[lapply(lapply(select(current_data,where(is.factor)),
                                                    unique),length) <2]
  if(length(single_value_factors) > 0) {
    current_data <- current_data |> select(-all_of(single_value_factors))
  }
  
  # Save dataset
  combo_name <- paste(gsub(" ", "_", as.character(current_type)),
                      gsub(" ", "_", as.character(current_diag)), sep = "_")
  
  data_file <- file.path(data_path, paste0(combo_name, ".RData"))
  save(current_data, file = data_file)
  
  # --- Summary for numeric variables ---
  num_vars <- current_data |> select(where(is.numeric))
  
  if (ncol(num_vars) > 0) {
    num_summary <- lapply(num_vars, summarize_numeric)
    num_summary <- do.call(rbind, num_summary)
    num_summary <- as.data.frame(num_summary)
    num_summary$Variable <- rownames(num_summary)
    rownames(num_summary) <- NULL
    num_summary <- num_summary |>
      select(Variable, everything())
  } else {
    num_summary <- NULL
  }
  numeric_file <- file.path(summaries_path, paste0("numeric_variables_", combo_name, ".csv")) 
  write.csv(num_summary, numeric_file, row.names = FALSE)
  
  
  # --- Summary for factor variables ---
  factor_vars <- current_data |> select(where(is.factor))
  
  if (ncol(factor_vars) > 0) {
    factor_summary <- lapply(names(factor_vars), function(v) {
      tbl <- table(factor_vars[[v]])
      out <- as.data.frame(t(tbl))
      out$Variable <- v
      out
    }) %>% bind_rows()
    
    # reorder Variable column first
    factor_summary <- factor_summary %>% select(Variable, everything())
    
  } else {
    factor_summary <- NULL
  }
  factor_file <- file.path(summaries_path, paste0("factor_variables_", combo_name, ".csv")) 
  write.csv(factor_summary, factor_file, row.names = FALSE)
  
  # Record counts
  counts_df <- rbind(counts_df, data.frame(
    type_diagnosis = combo_name,
    count = nrow(current_data),
    stringsAsFactors = FALSE
  ))
}

# Save counts and skipped log
counts_file <- file.path(results_path, "counts_by_type_diagnosis.csv")
write.csv(counts_df, counts_file, row.names = FALSE)

if(length(skipped_combinations) > 0){
  skipped_file <- file.path(results_path, "skipped_combinations.csv")
  write.csv(data.frame(skipped = skipped_combinations), skipped_file, row.names = FALSE)
}

cat("Data cleaning, splitting, and summaries completed.\n")
cat("Cleaned datasets saved to Data/, summaries to Results/Summaries/\n")
if(length(skipped_combinations) > 0) {
  cat("Skipped combinations recorded in:", skipped_file, "\n")
}
