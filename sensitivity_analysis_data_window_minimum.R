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

get_en_col <- function(sm) {
  if ("LogisticRegression_AUC" %in% names(sm)) "LogisticRegression_AUC"
  else "ElasticNet_AUC"
}

# Parse "[lo, hi]" strings into numeric, flag [NA,NA] / [0,1] / [1,1] as degenerate
# (mirrors step3_analysis-v2.R section 3b)
parse_ci <- function(x) suppressWarnings(as.numeric(strsplit(gsub("\\[|\\]", "", x), ",\\s*")[[1]]))
is_bad   <- function(ci) (is.na(ci[1]) && is.na(ci[2])) || isTRUE(all.equal(ci, c(0,1))) || isTRUE(all.equal(ci, c(1,1)))

# Load every analysis' participant summary for a given valid-window-minimum
# folder, null out any model whose own CI is degenerate, and return the
# StudyIDs to drop globally for that folder: participants with zero usable
# models in ANY config, or with only a single held-out observation in ANY
# config (mirrors step3_analysis-v2.R section 3b).
compute_drop_ids <- function(folder, analyses) {
  all_summaries <- list()
  for (an in analyses) {
    summary_path <- file.path(folder, an, "Summary",
                              paste0(an, "_participant_summary.csv"))
    if (file.exists(summary_path))
      all_summaries[[an]] <- read_csv(summary_path, show_col_types = FALSE)
  }
  if (length(all_summaries) == 0) return(character(0))

  all_summaries <- lapply(all_summaries, function(sm) {
    pairs <- list(c("RandomForest_AUC","RF_CI"),
                  c(get_en_col(sm), if ("LR_CI" %in% names(sm)) "LR_CI" else "EN_CI"),
                  c("GaussianProcess_AUC","GP_CI"),
                  c("Ensemble_AUC","Ens_CI"))
    for (p in pairs) {
      bad <- sapply(sm[[p[2]]], function(x) is_bad(parse_ci(x)))
      sm[[p[1]]][bad] <- NA
      sm[[p[2]]][bad] <- "[NA, NA]"
    }
    sm
  })

  auc_cols_by_sm <- lapply(all_summaries, function(sm) c("RandomForest_AUC", get_en_col(sm), "GaussianProcess_AUC", "Ensemble_AUC"))
  drop_ids_degenerate <- unique(unlist(lapply(names(all_summaries), function(an) {
    sm <- all_summaries[[an]]
    sm$StudyID[apply(sm[, auc_cols_by_sm[[an]]], 1, function(r) all(is.na(r)))]
  })))

  drop_ids_single_obs <- unique(unlist(lapply(names(all_summaries), function(an) {
    sm <- all_summaries[[an]]
    sm$StudyID[sm$NumObservations == 1]
  })))

  unique(c(drop_ids_degenerate, drop_ids_single_obs))
}

get_participant_auc <- function(folder, analysis_name, drop_ids = character(0)) {
  summary_path <- file.path(folder, analysis_name, "Summary",
                            paste0(analysis_name, "_participant_summary.csv"))
  if (!file.exists(summary_path)) return(NULL)

  summary_df <- read_csv(summary_path, show_col_types = FALSE)
  en_col <- get_en_col(summary_df)
  auc_cols <- c("RandomForest_AUC", en_col, "GaussianProcess_AUC", "Ensemble_AUC")

  # null out any model whose own CI is degenerate before taking the max
  pairs <- list(c("RandomForest_AUC","RF_CI"),
                c(en_col, if ("LR_CI" %in% names(summary_df)) "LR_CI" else "EN_CI"),
                c("GaussianProcess_AUC","GP_CI"),
                c("Ensemble_AUC","Ens_CI"))
  for (p in pairs) {
    bad <- sapply(summary_df[[p[2]]], function(x) is_bad(parse_ci(x)))
    summary_df[[p[1]]][bad] <- NA
  }

  summary_df %>%
    filter(!StudyID %in% drop_ids) %>%
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

# Drop IDs (degenerate CI / single-observation participants) computed once per
# valid-window-minimum folder, across all configurations within that folder.
drop_ids_by_folder <- setNames(
  lapply(min_folders, compute_drop_ids, analyses = analyses),
  min_folders
)

sensitivity_table <- bind_rows(lapply(analyses, function(analysis_name) {
  reference_auc <- get_participant_auc(reference_folder, analysis_name,
                                       drop_ids_by_folder[[reference_folder]])
  if (is.null(reference_auc)) return(NULL)

  bind_rows(lapply(min_folders, function(folder) {
    if (identical(folder, reference_folder)) return(NULL)

    comparison_auc <- get_participant_auc(folder, analysis_name,
                                          drop_ids_by_folder[[folder]])
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
drop_ids_auc_summary <- compute_drop_ids(auc_summary_folder, analyses)

auc_cutoff_summary <- bind_rows(lapply(auc_cutoffs, function(cutoff) {
  bind_rows(lapply(analyses, function(analysis_name) {
    auc_df <- get_participant_auc(auc_summary_folder, analysis_name, drop_ids_auc_summary)
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
