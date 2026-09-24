################################################################################
#
# Pain Increase Magnitude Analysis (standalone)
#
#   Rebuilds base_data from df_ema.csv using the same logic as
#   step2_imputation_and_modeling.R (Section 5a), restricts to the
#   participant IDs sourced from eligible_ids_local.R, and reports the
#   magnitude of pain increases (rows where pain_increasing == "Yes").
#
#   ASSUMPTIONS:
#     - eligible_ids_local.R lives in the same folder as this script and
#       defines a vector of StudyIDs under one of: eligible_ids, eligible,
#       eligible_studyids.
#     - df_ema.csv lives at the output_dir path below (same as
#       step1_create_datasets.R / step2_imputation_and_modeling.R produce).
#
################################################################################

library(readr)
library(dplyr)

## ---- 1. Paths -----------------------------------------------------------

source("paths_local.R")
ema_lookup_path <- file.path(output_dir, "df_ema.csv")

script_dir <- mice_reuse_path
eligible_ids_path <- file.path(script_dir, "eligible_ids_local.R")
if (!file.exists(eligible_ids_path)) {
  stop("eligible_ids_local.R not found at: ", eligible_ids_path,
       " — update script_dir/eligible_ids_path if it lives elsewhere.",
       call. = FALSE)
}

## ---- 2. Pull in eligible StudyIDs ---------------------------------------

source(eligible_ids_path)

eligible_vec <- if (exists("eligible_ids")) {
  eligible_ids
} else if (exists("eligible")) {
  eligible
} else if (exists("eligible_studyids")) {
  eligible_studyids
} else {
  stop("Could not find an eligible-IDs vector in eligible_ids_local.R. ",
       "Expected one of: eligible_ids, eligible, eligible_studyids.",
       call. = FALSE)
}
cat("Loaded", length(eligible_vec), "eligible StudyIDs from eligible_ids_local.R\n")

## ---- 3. Rebuild base_data (mirrors step2_imputation_and_modeling.R, Sec 5a) --

if (!file.exists(ema_lookup_path)) {
  stop("df_ema.csv not found at: ", ema_lookup_path,
       " — update ema_lookup_path at the top of this script.", call. = FALSE)
}

df_ema <- read_csv(ema_lookup_path, show_col_types = FALSE) %>%
  mutate(time_block = as.POSIXct(time_block, tz = "UTC"))

base_data <- df_ema %>%
  group_by(StudyID) %>%
  arrange(time_block) %>%
  mutate(
    days_since_first_ema = as.numeric(difftime(functional_date, min(functional_date), units = "days")),
    overall_pain_lag1 = lag(overall_pain, 1),
    diff              = overall_pain - overall_pain_lag1,
    pain_flag       = case_when(
      is.na(overall_pain_lag1) | is.na(diff) ~ NA_integer_,
      (diff > 0) ~ 1L, TRUE ~ 0L),
    total_transition  = max(cumsum(replace(pain_flag, is.na(pain_flag), 0))),
    n                 = n(),
    perc_transition   = total_transition / n,
    pain_increasing   = factor(
      case_when(
        is.na(pain_flag) ~ NA_character_,
        pain_flag == 1L  ~ "Yes",
        TRUE             ~ "No"
      ), levels = c("No","Yes")),
    overall_pain_lag2    = lag(overall_pain,    2),
    catastrophize_lag1   = lag(catastrophize,   1),
    catastrophize_lag2   = lag(catastrophize,   2),
    depress_lag1         = lag(depress,         1),
    depress_lag2         = lag(depress,         2),
    interference_lag1    = lag(interference,    1),
    interference_lag2    = lag(interference,    2),
    ema_missing_lag1     = lag(ema_missing,     1)
  ) %>%
  ungroup() %>%
  filter(!is.na(overall_pain)) %>%
  filter(perc_transition >= .2 & perc_transition <= .8)

## ---- 4. Restrict to eligible participants -------------------------------

base_data_eligible <- base_data %>% filter(StudyID %in% eligible_vec)

cat("base_data rows (all participants):     ", nrow(base_data), "\n")
cat("base_data rows (eligible participants):", nrow(base_data_eligible), "\n")
cat("Eligible IDs present in df_ema:         ",
    length(intersect(eligible_vec, unique(base_data_eligible$StudyID))),
    "/", length(eligible_vec), "\n")

## ---- 5. Prevalence: how many/what proportion of rows are flagged --------

pain_flag_summary <- base_data_eligible %>%
  filter(!is.na(pain_increasing)) %>%
  summarise(
    n_total      = n(),
    n_increase   = sum(pain_increasing == "Yes"),
    pct_increase = round(100 * n_increase / n_total, 2)
  )
cat("\n--- Pain increase prevalence (eligible participants) ---\n")
print(pain_flag_summary)

## ---- 6. Magnitude: size of the increase itself ---------------------------

increase_magnitude <- base_data_eligible %>%
  filter(pain_increasing == "Yes") %>%
  summarise(
    n           = n(),
    mean_diff   = mean(diff, na.rm = TRUE),
    median_diff = median(diff, na.rm = TRUE),
    sd_diff     = sd(diff, na.rm = TRUE),
    min_diff    = min(diff, na.rm = TRUE),
    max_diff    = max(diff, na.rm = TRUE),
    q25         = quantile(diff, 0.25, na.rm = TRUE),
    q75         = quantile(diff, 0.75, na.rm = TRUE)
  )
cat("\n--- Magnitude of pain increases (eligible participants) ---\n")
print(increase_magnitude)
# 
# ## ---- 7. Per-participant breakdown ----------------------------------------
# 
# increase_magnitude_by_pid <- base_data_eligible %>%
#   filter(pain_increasing == "Yes") %>%
#   group_by(StudyID) %>%
#   summarise(
#     n_increases = n(),
#     mean_diff   = mean(diff, na.rm = TRUE),
#     median_diff = median(diff, na.rm = TRUE),
#     max_diff    = max(diff, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(mean_diff))
# cat("\n--- Per-participant increase magnitude (eligible participants) ---\n")
# print(increase_magnitude_by_pid)

## ---- 8. Save summary + histogram ------------------------------------------

# out_csv <- file.path(script_dir, "pain_increase_magnitude_by_pid.csv")
# write_csv(increase_magnitude_by_pid, out_csv)
# cat("\nSaved per-participant summary to:", out_csv, "\n")
# 
# out_png <- file.path(script_dir, "pain_increase_magnitude_hist.png")
# png(out_png, width = 800, height = 600)
# hist(base_data_eligible$diff[base_data_eligible$pain_increasing == "Yes"],
#      main = "Magnitude of pain increases (eligible participants)",
#      xlab = "Increase in overall_pain (diff = overall_pain - lag1)", breaks = 20)
# dev.off()
# cat("Saved histogram to:", out_png, "\n")