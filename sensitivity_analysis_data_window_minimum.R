################################################################################
#
# Supplemental Analysis: Valid-Window Minimum Sensitivity
#
#   Aggregates pooled performance tables across valid-window minimums,
#   compares participant-level best AUC values against a reference minimum,
#   and summarizes the number of participants exceeding selected AUC cutoffs.
#
################################################################################

library(readr)
library(dplyr)


################################################################################
# 0 - CONFIGURATION
################################################################################

pipeline_dir <- "/Users/f0085f6/Desktop/Frumkin_lab/personalizedPainPrediction_paper/PPP_Project/Output/pipeline/simple_pain_outcome"
reference_min <- "74ValidWindows"
auc_summary_min <- "26ValidWindows"

analyses <- c("A1_FitbitOnly", "A2_LagEMAOnly", "A3_Combined", "A4_Lag1PainOnly")
auc_cutoffs <- c(0.7, 0.8)


################################################################################
# 1 - HELPER FUNCTIONS
################################################################################

list_min_folders <- function(pipeline_dir) {
  folders <- list.dirs(pipeline_dir, full.names = TRUE, recursive = FALSE)
  folders[grepl("ValidWindows", basename(folders))]
}

get_participant_auc <- function(folder, analysis_name) {
  summary_path <- file.path(folder, analysis_name, "Summary",
                            paste0(analysis_name, "_participant_summary.csv"))
  if (!file.exists(summary_path)) return(NULL)

  summary_df <- read_csv(summary_path, show_col_types = FALSE)
  en_col <- if ("LogisticRegression_AUC" %in% names(summary_df)) {
    "LogisticRegression_AUC"
  } else {
    "ElasticNet_AUC"
  }
  auc_cols <- c("RandomForest_AUC", en_col, "GaussianProcess_AUC", "Ensemble_AUC")

  summary_df %>%
    mutate(Best_AUC = apply(select(., all_of(auc_cols)), 1, max, na.rm = TRUE)) %>%
    select(StudyID, Best_AUC)
}

format_pvalue <- function(p, digits = 3) {
  ifelse(p < 0.001, "<.001", formatC(round(p, digits), format = "f", digits = digits))
}


################################################################################
# 2 - COMBINE TABLE 2 ACROSS MINIMUMS
################################################################################

min_folders <- list_min_folders(pipeline_dir)

cat("Found minimum folders:\n")
print(basename(min_folders))

combined_table2 <- bind_rows(lapply(min_folders, function(folder) {
  table_path <- file.path(folder, "Table2_performance_metrics_pooled.csv")
  if (!file.exists(table_path)) {
    cat("  Missing:", basename(folder), "\n")
    return(NULL)
  }

  read_csv(table_path, show_col_types = FALSE) %>%
    mutate(Minimum = basename(folder))
})) %>%
  select(-any_of("Note")) %>%
  select(`Feature Configuration`, Minimum, everything()) %>%
  arrange(`Feature Configuration`, Minimum)

write_csv(combined_table2, file.path(pipeline_dir, "Table2_combined_minimums.csv"))
cat("\nSaved: Table2_combined_minimums.csv\n")


################################################################################
# 3 - AUC SENSITIVITY AGAINST REFERENCE MINIMUM
################################################################################

reference_folder <- min_folders[basename(min_folders) == reference_min]

sensitivity_table <- bind_rows(lapply(analyses, function(analysis_name) {
  reference_auc <- get_participant_auc(reference_folder, analysis_name)
  if (is.null(reference_auc)) return(NULL)

  bind_rows(lapply(min_folders, function(folder) {
    if (identical(folder, reference_folder)) return(NULL)

    comparison_auc <- get_participant_auc(folder, analysis_name)
    if (is.null(comparison_auc)) return(NULL)

    tt <- t.test(reference_auc$Best_AUC, comparison_auc$Best_AUC)

    data.frame(
      `Feature Configuration`  = analysis_name,
      `Valid Window Minimum` = as.integer(gsub("ValidWindows", "", basename(folder))),
      `N (Reference)`          = nrow(reference_auc),
      `N (Other)`              = nrow(comparison_auc),
      `Mean AUC (Reference)`   = round(mean(reference_auc$Best_AUC), 3),
      `Mean AUC (Other)`       = round(mean(comparison_auc$Best_AUC), 3),
      Difference               = round(mean(reference_auc$Best_AUC) - mean(comparison_auc$Best_AUC), 3),
      `Unadjusted p-value`     = tt$p.value,
      check.names = FALSE
    )
  }))
})) %>%
  group_by(`Feature Configuration`) %>%
  mutate(`Adjusted p-value` = p.adjust(`Unadjusted p-value`, method = "BH")) %>%
  ungroup() %>%
  mutate(
    `Unadjusted p-value` = format_pvalue(`Unadjusted p-value`),
    `Adjusted p-value` = format_pvalue(`Adjusted p-value`)
  )

write_csv(sensitivity_table, file.path(pipeline_dir, "AUC_minimum_sensitivity.csv"))


################################################################################
# 4 - PARTICIPANT COUNTS ABOVE AUC CUTOFFS
################################################################################

auc_summary_folder <- file.path(pipeline_dir, auc_summary_min)

auc_cutoff_summary <- bind_rows(lapply(auc_cutoffs, function(cutoff) {
  bind_rows(lapply(analyses, function(analysis_name) {
    auc_df <- get_participant_auc(auc_summary_folder, analysis_name)
    if (is.null(auc_df)) return(NULL)

    data.frame(
      `Feature Configuration` = analysis_name,
      AUC_Cutoff = cutoff,
      N_total = nrow(auc_df),
      N_above_cutoff = sum(auc_df$Best_AUC >= cutoff, na.rm = TRUE),
      Pct_above_cutoff = round(mean(auc_df$Best_AUC >= cutoff, na.rm = TRUE) * 100, 1),
      check.names = FALSE
    )
  }))
}))

print(auc_cutoff_summary)
