################################################################################
#
# Supplemental Tables
#
#   Produces participant demographics and predictor missingness tables for the
#   model-eligible sample.
#
################################################################################

library(readr)
library(dplyr)
library(lubridate)
library(data.table)
library(tibble)
library(tableone)


################################################################################
# 0 - CONFIGURATION
################################################################################

source("paths_local.R")

demographics_path <- file.path(data_dir, "redcap_demographics_all.csv")
ema_source_path <- file.path(data_dir, "allEMAdata-v5.csv")
expanded_path <- file.path(output_dir, "df_expanded.csv")
ema_lookup_path <- file.path(output_dir, "df_ema.csv")
eligible_ids_path <- "eligible_ids_local.R"

source(eligible_ids_path)

demographics_vars <- c(
  "age",
  "gender",
  "race",
  "ethnicity",
  "time0_pain_intesity_score",
  "time0_pain_inter_score",
  "time0_phys_func_score"
)

demographics_factor_vars <- c("gender", "race", "ethnicity")


################################################################################
# 1 - PARTICIPANT DEMOGRAPHICS
################################################################################

# Recode categorical fields to readable labels before tabulation.
df_demographics <- read_csv(demographics_path, show_col_types = FALSE) %>%
  mutate(
    gender = factor(gender, levels = c(1, 2, 3),
                    labels = c("Male", "Female", "Other")),
    race = factor(race, levels = c(1, 2, 3, 4, 5),
                  labels = c("White", "African American/Black", "Asian",
                             "Native American", "Other")),
    ethnicity = factor(ethnicity, levels = c(1, 2),
                       labels = c("Hispanic/Latino", "Not Hispanic/Latino"))
  ) %>%
  filter(StudyID %in% eligible_ids) %>%
  select(StudyID, all_of(demographics_vars))

demographics_table <- CreateTableOne(
  vars = demographics_vars,
  factorVars = demographics_factor_vars,
  data = df_demographics
)

print(demographics_table, showAllLevels = TRUE, formatOptions = list(big.mark = ","))


################################################################################
# 2 - PAIN OUTCOME DISTRIBUTION & PREDICTOR MISSINGNESS
################################################################################

# Use the pre-modeling EMA file to summarize missingness before imputation.
df_ema_raw <- read_csv(ema_source_path, show_col_types = FALSE) %>%
  arrange(StudyID, Response.Time.New) %>%
  filter(StudyID %in% eligible_ids, !is.na(Notification.Time)) %>%
  mutate(
    pain_missing = as.integer(is.na(pain_intensity)),
    cat_missing = as.integer(is.na(catastrophizing)),
    interference_missing = as.integer(is.na(pain_interference)),
    depression_missing = as.integer(is.na(depressed)),
    opioid_missing = as.integer(is.na(opioid))
  )

ema_missingness_table <- tibble(
  Variable = c("Pain intensity", "Catastrophizing", "Pain interference", "Depression"),
  N_Missing = c(
    sum(df_ema_raw$pain_missing),
    sum(df_ema_raw$cat_missing),
    sum(df_ema_raw$interference_missing),
    sum(df_ema_raw$depression_missing)
  ),
  Pct_Missing = round(c(
    mean(df_ema_raw$pain_missing),
    mean(df_ema_raw$cat_missing),
    mean(df_ema_raw$interference_missing),
    mean(df_ema_raw$depression_missing)
  ) * 100, 2)
)

opioid_missingness <- df_ema_raw %>%
  group_by(StudyID, functional_date) %>%
  summarise(opioid_missing_day = unique(opioid_missing), .groups = "drop") %>%
  summarise(
    Variable = "Opioid use",
    N_Missing = sum(opioid_missing_day),
    Pct_Missing = round(mean(opioid_missing_day) * 100, 2)
  )

df_expanded <- read_csv(expanded_path, show_col_types = FALSE) %>%
  mutate(time_block = as.POSIXct(time_block, tz = "UTC"))

df_ema <- read_csv(ema_lookup_path, show_col_types = FALSE) %>%
  mutate(time_block = as.POSIXct(time_block, tz = "UTC"))

# Applies the same eligibility filter used in the modeling pipeline.
base_data <- df_ema %>%
  group_by(StudyID) %>%
  arrange(time_block) %>%
  mutate(
    overall_pain_lag1 = lag(overall_pain, 1),
    diff = overall_pain - overall_pain_lag1,
    pain_flag = case_when(
      is.na(overall_pain_lag1) | is.na(diff) ~ NA_integer_,
      diff > 0 ~ 1L,
      TRUE ~ 0L
    ),
    total_transition = max(cumsum(replace(pain_flag, is.na(pain_flag), 0))),
    n = n(),
    perc_transition = total_transition / n,
    pain_increasing = factor(
      case_when(
        is.na(pain_flag) ~ NA_character_,
        pain_flag == 1L ~ "Yes",
        TRUE ~ "No"
      ),
      levels = c("No", "Yes")
    )
  ) %>%
  ungroup() %>%
  filter(perc_transition >= 0.2 & perc_transition <= 0.8)

ema_windows_full <- base_data %>%
  filter(
    StudyID %in% eligible_ids,
    !is.na(pain_increasing),
    !is.na(catastrophize),
    !is.na(depress),
    !is.na(interference)
  )

# Get classification distribution table
avg_diff_by_class <- ema_windows_full %>%
  group_by(pain_increasing) %>%
  summarise(Avg_Diff = round(mean(diff, na.rm = TRUE), 2), .groups = "drop")

pain_increasing_table <- tibble(
  Class = c("No", "Yes"),
  N = as.integer(table(ema_windows_full$pain_increasing)),
  Pct = round(as.numeric(prop.table(table(ema_windows_full$pain_increasing))) * 100, 1)
) %>%
  left_join(avg_diff_by_class, by = c("Class" = "pain_increasing")) %>%
  bind_rows(tibble(
    Class = "Total",
    N = sum(.$N),
    Pct = sum(.$Pct),
    Avg_Diff = round(mean(ema_windows_full$diff, na.rm = TRUE), 2)
  ))

cat("\nPain-increasing class distribution (eligible sample):\n")
print(pain_increasing_table)

ema_windows <- ema_windows_full %>%
  select(StudyID, time_block, functional_date) %>%
  mutate(
    window_start = time_block - hours(1),
    window_end = time_block
  )

sensor_dt <- as.data.table(df_expanded)[
  ,
  .(StudyID, time_block, hr_missing, steps_missing, sleep_missing)
]
window_dt <- as.data.table(ema_windows)

# Compute 1-hour pre-EMA sensor-window missingness.
matched <- sensor_dt[
  window_dt,
  on = .(StudyID, time_block > window_start, time_block <= window_end),
  .(
    StudyID,
    ema_time = i.time_block,
    functional_date = i.functional_date,
    hr_missing,
    steps_missing,
    sleep_missing
  ),
  nomatch = 0
]

window_missingness <- matched[
  ,
  .(
    hr_missing_win = mean(hr_missing, na.rm = TRUE),
    steps_missing_win = mean(steps_missing, na.rm = TRUE),
    sleep_missing_win = as.integer(mean(sleep_missing, na.rm = TRUE) > 0)
  ),
  by = .(StudyID, ema_time, functional_date)
]

all_windows <- window_dt[, .(StudyID, ema_time = time_block, functional_date)]
window_missingness <- merge(
  all_windows,
  window_missingness,
  by = c("StudyID", "ema_time", "functional_date"),
  all.x = TRUE
)

window_missingness[is.na(hr_missing_win), hr_missing_win := 1]
window_missingness[is.na(steps_missing_win), steps_missing_win := 1]
window_missingness[is.na(sleep_missing_win), sleep_missing_win := 1L]

sensor_missingness_table <- tibble(
  Variable = c("Heart rate (window)", "Steps (window)", "Sleep (window)"),
  N_Missing = c(
    sum(window_missingness$hr_missing_win),
    sum(window_missingness$steps_missing_win),
    sum(window_missingness$sleep_missing_win)
  ),
  Pct_Missing = round(c(
    mean(window_missingness$hr_missing_win),
    mean(window_missingness$steps_missing_win),
    mean(window_missingness$sleep_missing_win)
  ) * 100, 2)
)

predictor_missingness_table <- bind_rows(
  ema_missingness_table,
  opioid_missingness,
  sensor_missingness_table
)

print(predictor_missingness_table)


################################################################################
# 3 - REPLICATE STEP 2'S BASE_DATA + DATA WINDOWS PER PARTICIPANT
################################################################################

# Adds the remaining EMA lag columns to base_data.
base_data_full <- base_data %>%
  group_by(StudyID) %>%
  arrange(time_block) %>%
  mutate(
    days_since_first_ema = as.numeric(difftime(functional_date, min(functional_date), units = "days")),
    overall_pain_lag2   = lag(overall_pain,  2),
    catastrophize_lag1  = lag(catastrophize, 1),
    catastrophize_lag2  = lag(catastrophize, 2),
    depress_lag1        = lag(depress,       1),
    depress_lag2        = lag(depress,       2),
    interference_lag1   = lag(interference,  1),
    interference_lag2   = lag(interference,  2),
    ema_missing_lag1     = lag(ema_missing,   1)
  ) %>%
  ungroup() %>%
  filter(!is.na(overall_pain))

# Rows of data ("data windows") per eligible participant
ema_complete_all <- base_data_full %>%
  filter(
    StudyID %in% eligible_ids,
    !is.na(pain_increasing),
    !is.na(overall_pain_lag2),
    !is.na(catastrophize_lag1), !is.na(catastrophize_lag2),
    !is.na(depress_lag1),       !is.na(depress_lag2),
    !is.na(interference_lag1),  !is.na(interference_lag2),
    !is.na(ema_missing_lag1)
  )

windows_per_participant <- ema_complete_all %>%
  count(StudyID, name = "N_Windows")

cat("\nData windows (rows) per participant (N =", nrow(windows_per_participant), "):\n")
cat("  Mean:  ", round(mean(windows_per_participant$N_Windows), 1), "\n")
cat("  SD:    ", round(sd(windows_per_participant$N_Windows), 1), "\n")
cat("  Median:", median(windows_per_participant$N_Windows), "\n")
cat("  Range: ", min(windows_per_participant$N_Windows), "-",
    max(windows_per_participant$N_Windows), "\n")

print(windows_per_participant %>% arrange(N_Windows))
