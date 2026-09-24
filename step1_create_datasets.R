################################################################################
#
# Personalized Pain Prediction (PPP) — Step 1: Dataset Creation
#
#   PART 1 (Sections 1-4):  Merge raw EMA, Fitbit, and sleep data; flag
#                           missingness; add time descriptors; aggregate
#                           into 5-minute time blocks.
#                           -> writes allMergedData_5minBlocks-v4.csv
#
#   PART 2 (Sections 5-10): Stage 1 preprocessing for modeling. Loads the
#                           5-minute block data from Part 1, expands each
#                           participant to a full 5-minute time grid, and
#                           builds an EMA lookup table merged with weather
#                           data.
#                           -> writes df_expanded.csv and df_ema.csv
#
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(lubridate)


################################################################################
# 0 - CONFIGURATION
################################################################################

source("paths_local.R")
report_dir <- file.path(output_dir, "reports")

for (d in c(output_dir, report_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

fitbit_path  <- file.path(data_dir, "allFitbitData-v3.csv")
sleep_path   <- file.path(data_dir, "allSleepData-v2.csv")
ema_raw_path <- file.path(data_dir, "Raw EMA Data", "EMA_data_long-v2.csv")
demo_path    <- file.path(data_dir, "Old datasets", "merged_data.csv")

weather_path <- file.path(data_dir, "weather_chunks",
                          "st louis, MO, USA 2021-02-01 to 2023-08-01.csv")

ema_output_path    <- file.path(data_dir, "allEMAdata-v5.csv")
merged_output_path <- file.path(data_dir, "allMergedData-v5.csv")
blocks_path        <- file.path(data_dir, "allMergedData_5minBlocks-v4.csv")

expanded_path   <- file.path(output_dir, "df_expanded.csv")
ema_lookup_path <- file.path(output_dir, "df_ema.csv")

# EMA observations before this hour are attributed to the previous study day.
study_day_cutoff_hour <- 3


################################################################################
# PART 1 - BUILD MERGED 5-MINUTE-BLOCK DATASET
################################################################################

###############################################################################
# 1 - MERGE EMA, FITBIT, AND SLEEP DATA
###############################################################################

cat("=== Part 1: Building merged dataset ===\n")
cat("\n[1/4] Loading and merging EMA, Fitbit, and sleep data...\n")

# Fitbit heart rate values of -1 indicate missing sensor readings.
df_fitbit <- read_csv(fitbit_path, show_col_types = FALSE) %>%
  separate(timestamp, c("date", "time"), sep = " ", remove = FALSE) %>%
  mutate(
    hr        = ifelse(hr == -1, NA, hr),
    steps     = ifelse(hr == -1, NA, steps),
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    date      = as.Date(date),
    StudyID   = str_remove(StudyID, "_v1")
  ) %>%
  group_by(StudyID, timestamp) %>%
  arrange(hr, steps) %>%
  slice(1) %>%
  ungroup()

cat("  Fitbit rows:", nrow(df_fitbit), "\n")

df_sleep_all <- read_csv(sleep_path, show_col_types = FALSE) %>%
  mutate(
    date    = as.Date(date),
    StudyID = str_remove(StudyID, "_v1")
  )

df_sleep <- df_sleep_all %>%
  select(-c(main_startTime, time_series, totalSleepRecords,
            main_restlessCounts, main_remCounts))

cat("  Sleep rows:", nrow(df_sleep), "\n")

df_sensor_merged <- full_join(df_fitbit, df_sleep, by = c("StudyID", "date"))

cat("  Combined sensor rows:", nrow(df_sensor_merged), "\n")

df_ema_source <- read_csv(ema_raw_path, show_col_types = FALSE)

df_ema_processed <- df_ema_source %>%
  filter(str_detect(composite_label, "preop")) %>%
  mutate(composite_label = str_remove(composite_label, "preop.")) %>%
  mutate(value = case_when(
    composite_label == "pain_intensity" ~ composite_max,
    TRUE                                ~ composite_mean
  )) %>%
  select(-c(composite_max, composite_mean)) %>%
  pivot_wider(names_from = composite_label, values_from = value) %>%
  select(-`NA`) %>%
  mutate(
    Response.Time.New = as.POSIXct(Response.Time.New, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    Notification.Time  = as.POSIXct(Notification.Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ) %>%
  mutate(ts_for_date = if_else(is.na(Response.Time.New), Notification.Time, Response.Time.New)) %>%
  separate(ts_for_date, c("date", "time"), sep = " ", remove = FALSE) %>%
  separate(time, c("h"), extra = "drop") %>%
  group_by(StudyID) %>%
  arrange(StudyID, Notification.Time, Response.Time.New) %>%
  mutate(
    h    = as.numeric(h),
    h    = ifelse(h < 3, h + 24, h),
    date = as.Date(date),
    day  = as.numeric(date - first(date)),
    time = dense_rank(h),
    time_s = as.numeric(Response.Time.New - first(na.omit(Response.Time.New))),
    functional_date = as.Date(Notification.Time - hours(study_day_cutoff_hour))
  ) %>%
  select(-c(h, ts_for_date)) %>%
  filter(!is.na(Notification.Time)) %>%
  arrange(StudyID, functional_date, time) %>%
  group_by(StudyID, functional_date) %>%
  fill(opioid, .direction = "up") %>%
  ungroup()

# SPINE2166 contains duplicated EMA schedules on some days.
df_ema_clean <- df_ema_processed %>%
  arrange(StudyID, functional_date, time) %>%
  group_by(StudyID, functional_date) %>%
  filter(StudyID != "SPINE2166" | n() <= 5 | row_number() %% 2 == 1) %>%
  ungroup()

cat("  EMA rows (post-cleaning):", nrow(df_ema_clean), "\n")

df_ema_clean %>%
  mutate(Response.Time.New = format(Response.Time.New, "%Y-%m-%d %H:%M:%S")) %>%
  write_csv(ema_output_path)

df_ema_sensor <- full_join(
  df_ema_clean, df_sensor_merged,
  by = c("StudyID", "Response.Time.New" = "timestamp", "date")
)

df_demo <- read_csv(demo_path, show_col_types = FALSE) %>%
  select(c(StudyID, Surgery_date, Age, Sex)) %>%
  distinct()

df_merged <- df_ema_sensor %>%
  left_join(df_demo, by = "StudyID") %>%
  mutate(Surgery_date = mdy(Surgery_date)) %>%
  filter(date < Surgery_date | is.na(Surgery_date))

df_merged %>%
  mutate(Response.Time.New = format(Response.Time.New, "%Y-%m-%d %H:%M:%S")) %>%
  write_csv(merged_output_path)

cat("  Merged rows:", nrow(df_merged), "\n")


###############################################################################
# 2 - FLAG MISSING VARIABLES
###############################################################################

cat("\n[2/4] Flagging missing variables...\n")

df_merged <- read_csv(merged_output_path, show_col_types = FALSE) %>%
  arrange(StudyID, Response.Time.New)

df_missing_flagged <- df_merged %>%
  mutate(
    ema_missing   = ifelse(!is.na(Notification.Time) & is.na(pain_intensity), 1, 0),
    hr_missing    = ifelse(is.na(hr), 1, 0),
    steps_missing = ifelse(is.na(steps), 1, 0),
    sleep_missing = ifelse(is.na(totalMinutesAsleep), 1, 0)
  )

cat("  ema_missing:",   sum(df_missing_flagged$ema_missing),   "\n")
cat("  hr_missing:",    sum(df_missing_flagged$hr_missing),    "\n")
cat("  steps_missing:", sum(df_missing_flagged$steps_missing), "\n")
cat("  sleep_missing:", sum(df_missing_flagged$sleep_missing), "\n")


###############################################################################
# 3 - ADD TIME DESCRIPTORS
###############################################################################

cat("\n[3/4] Adding time descriptors...\n")

df_time <- df_missing_flagged %>%
  group_by(StudyID) %>%
  mutate(
    hour = hour(Response.Time.New),
    time_of_day = case_when(
      hour >= 5  & hour < 12 ~ "Morning",
      hour >= 12 & hour < 17 ~ "Afternoon",
      hour >= 17 & hour < 21 ~ "Evening",
      hour >= 21 | hour < 5  ~ "Night",
      TRUE                   ~ NA_character_
    ),
    day_of_week = weekdays(date),
    is_weekend  = ifelse(day_of_week %in% c("Saturday", "Sunday"), TRUE, FALSE)
  ) %>%
  ungroup()


###############################################################################
# 4 - GROUP DATA INTO 5-MINUTE BLOCKS
###############################################################################

cat("\n[4/4] Grouping into 5-minute blocks...\n")

# Steps >= 200 in a single 5-minute interval are treated as sensor artifacts.
df_blocks <- df_time %>%
  filter(!is.na(Response.Time.New)) %>%
  filter(steps < 200 | is.na(steps)) %>%
  mutate(time_block = floor_date(Response.Time.New, unit = "5 minutes")) %>%
  group_by(StudyID, date, Surgery_date, Age, Sex, hour, time_of_day,
           day_of_week, is_weekend, time_block) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    total_steps = sum(steps, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(ema_missing = ifelse(ema_missing > 0, 1, 0)) %>%
  ungroup()

df_blocks[sapply(df_blocks, is.infinite)] <- NA

df_blocks %>%
  mutate(time_block = format(time_block, "%Y-%m-%d %H:%M:%S")) %>%
  write_csv(blocks_path)

cat("  5-minute block rows:", nrow(df_blocks), "\n")
cat("\n=== Part 1 complete ===\n")
cat("Wrote:", blocks_path, "\n")


################################################################################
# PART 2 - STAGE 1 PREPROCESSING
################################################################################

cat("\n=== Part 2 / Stage 1: Preprocessing ===\n")


###############################################################################
# 5 - LOAD RAW DATA & RENAME EMA VARIABLES
###############################################################################

cat("\n[1/5] Loading raw 5-minute block data...\n")

df_raw <- read_csv(blocks_path, show_col_types = FALSE) %>%
  mutate(
    time_block = as.POSIXct(time_block, tz = "UTC"),
    date       = as.Date(time_block)
  ) %>%
  arrange(StudyID, time_block)

cat("  Raw rows:", nrow(df_raw),
    "| Participants:", length(unique(df_raw$StudyID)), "\n")

# -- Rename EMA variables to pipeline-standard names --------------------------
# Raw name          -> Pipeline name
# pain_intensity    -> overall_pain
# catastrophizing   -> catastrophize
# depressed         -> depress
# pain_interference -> interference
# opioid            -> opioid_num
df_raw <- df_raw %>%
  rename(
    overall_pain  = pain_intensity,
    catastrophize = catastrophizing,
    depress       = depressed,
    interference  = pain_interference,
    opioid_num    = opioid
  )

cat("  Columns renamed to pipeline-standard names.\n")
cat("  EMA variables present:",
    paste(intersect(c("overall_pain", "catastrophize", "depress",
                      "interference", "opioid_num"), names(df_raw)),
          collapse = ", "), "\n")


###############################################################################
# 6 - ADD TIME DESCRIPTORS (GRID-CONSISTENT)
#
# Missingness flags (hr_missing, steps_missing, sleep_missing, ema_missing)
# already exist in the raw file as binary 0/1 - not recomputed here.
###############################################################################

cat("\n[2/5] Adding time descriptors...\n")

df_time_desc <- df_raw %>%
  mutate(
    # Recompute hour/minute from time_block.
    hour            = hour(time_block),
    minute          = minute(time_block),
    # functional_date shifts time_block back by the cutoff hour so pre-cutoff
    # rows belong to the previous study day.
    functional_date = as.Date(time_block - hours(study_day_cutoff_hour)),
    time_of_day     = case_when(
      hour >= 5  & hour < 12 ~ "Morning",
      hour >= 12 & hour < 17 ~ "Afternoon",
      hour >= 17 & hour < 21 ~ "Evening",
      TRUE                   ~ "Night"
    ),
    day_of_week = weekdays(functional_date),
    is_weekend  = day_of_week %in% c("Saturday", "Sunday")
  )


###############################################################################
# 7 - IDENTIFY NO-SENSOR PARTICIPANTS
###############################################################################

no_sensor_ids <- df_time_desc %>%
  group_by(StudyID) %>%
  summarise(
    no_hr    = all(is.na(hr)),
    no_steps = all(is.na(steps)),
    no_sleep = all(is.na(totalMinutesAsleep)),
    .groups  = "drop"
  ) %>%
  filter(no_hr & no_steps & no_sleep) %>%
  pull(StudyID)

cat("  No-sensor participants (excluded from imputation):",
    length(no_sensor_ids), "\n")
if (length(no_sensor_ids) > 0) print(no_sensor_ids)

df_sensor    <- df_time_desc %>% filter(!StudyID %in% no_sensor_ids)
df_no_sensor <- df_time_desc %>% filter( StudyID %in% no_sensor_ids)

# NOTE: no-sensor participants are excluded from grid expansion and therefore
# from df_ema.
cat("  Sensor-eligible participants for grid expansion:",
    length(unique(df_sensor$StudyID)), "\n")


###############################################################################
# 8 - EXPAND TO FULL 5-MIN GRID
#
# Grid spans first to last EMA timestamp per participant.
# An EMA row is any row where a notification was sent:
#   !is.na(overall_pain) | ema_missing > 0
###############################################################################

cat("\n[3/5] Expanding to full 5-min grid (this may take a few minutes)...\n")

ema_range <- df_sensor %>%
  filter(!is.na(overall_pain) | ema_missing > 0) %>%
  group_by(StudyID) %>%
  summarise(
    first_ema = min(time_block, na.rm = TRUE),
    last_ema  = max(time_block, na.rm = TRUE),
    .groups   = "drop"
  )

cat("  Participants with at least one EMA observation:", nrow(ema_range), "\n")

full_grid <- ema_range %>%
  rowwise() %>%
  mutate(time_block = list(seq(first_ema, last_ema, by = "5 mins"))) %>%
  unnest(time_block) %>%
  select(StudyID, time_block)

df_expanded <- full_grid %>%
  left_join(df_sensor, by = c("StudyID", "time_block")) %>%
  group_by(StudyID) %>%
  fill(Surgery_date, Age, Sex, .direction = "downup") %>%
  mutate(
    date            = as.Date(time_block),
    hour            = hour(time_block),
    minute          = minute(time_block),
    # functional_date shifts time_block back by the cutoff hour so pre-cutoff
    # rows belong to the previous study day.
    functional_date = as.Date(time_block - hours(study_day_cutoff_hour)),
    day             = as.integer(functional_date - min(functional_date, na.rm = TRUE)),
    time_of_day     = case_when(
      hour >= 5  & hour < 12 ~ "Morning",
      hour >= 12 & hour < 17 ~ "Afternoon",
      hour >= 17 & hour < 21 ~ "Evening",
      TRUE                   ~ "Night"
    ),
    day_of_week = weekdays(functional_date),
    is_weekend  = day_of_week %in% c("Saturday", "Sunday")
  ) %>%
  # Broadcasts the observed sleep value across all rows in the same
  # participant study-day.
  group_by(StudyID, functional_date) %>%
  mutate(
    totalMinutesAsleep = ifelse(
      is.na(totalMinutesAsleep),
      first(na.omit(totalMinutesAsleep)),
      totalMinutesAsleep)
  ) %>%
  ungroup() %>%
  mutate(
    # ema_missing: grid-expansion rows have NA (not a notification row) -> 0
    ema_missing   = ifelse(is.na(ema_missing),   0L, ema_missing),
    hr_missing    = ifelse(is.na(hr_missing),    1L, hr_missing),
    steps_missing = ifelse(is.na(steps_missing), 1L, steps_missing),
    sleep_missing = ifelse(is.na(sleep_missing),
                           ifelse(is.na(totalMinutesAsleep), 1L, 0L),
                           sleep_missing)
  )

cat("  Grid-expanded rows:", nrow(df_expanded), "\n")

df_expanded %>%
  mutate(time_block = format(time_block, "%Y-%m-%d %H:%M:%S")) %>%
  write_csv(expanded_path)

cat("  Saved:", expanded_path, "\n")


###############################################################################
# 9 - EMA LOOKUP: WEATHER
###############################################################################

cat("\n[4/5] Extracting EMA observations...\n")

# -- 9a: one row per actual EMA observation --------------------------------------
# EMA rows: notification was sent (response recorded OR missed)
df_ema_obs <- df_expanded %>%
  filter(!is.na(overall_pain) | ema_missing > 0) %>%
  mutate(
    Response.Date = as.Date(time_block),
    Response.Hour = hour(time_block)
  ) %>%
  distinct(StudyID, time_block, .keep_all = TRUE) %>%
  arrange(StudyID, Response.Date, Response.Hour)

cat("  EMA observations:", nrow(df_ema_obs), "\n")
cat("  With response (!is.na(overall_pain)):",
    sum(!is.na(df_ema_obs$overall_pain)), "\n")
cat("  Missed (ema_missing == 1, no response):",
    sum(is.na(df_ema_obs$overall_pain) & df_ema_obs$ema_missing == 1), "\n")


# -- 9b: weather merge ------------------------------------------------------------
cat("\n[5/5] Merging weather data...\n")

df_weather <- read_csv(weather_path, show_col_types = FALSE) %>%
  select(datetime, temp, humidity, precip, sealevelpressure) %>%
  rename(
    weather_temp     = temp,
    weather_precip   = precip,
    weather_humidity = humidity,
    weather_pressure = sealevelpressure
  ) %>%
  mutate(date = as.Date(datetime)) %>%
  distinct(date, .keep_all = TRUE) %>%
  select(date, weather_temp, weather_humidity, weather_pressure, weather_precip)

df_ema_weather <- df_ema_obs %>%
  mutate(
    date            = as.Date(time_block),
    functional_date = as.Date(time_block - hours(study_day_cutoff_hour))
  ) %>%
  left_join(df_weather, by = c("functional_date" = "date"))

weather_miss <- tibble(
  Variable    = c("weather_temp", "weather_humidity",
                  "weather_pressure", "weather_precip"),
  Missing_Pct = sapply(c("weather_temp", "weather_humidity",
                         "weather_pressure", "weather_precip"),
                       function(v) mean(is.na(df_ema_weather[[v]])) * 100)
)
cat("  Weather missingness (%):\n"); print(weather_miss)


# Keep only EMA-level variables here; sensor window features are rebuilt in
# step2_imputation_and_modeling.R.
ema_cols <- c(
  "StudyID", "time_block",
  "overall_pain", "catastrophize", "depress", "interference", "opioid_num",
  "weather_temp", "weather_humidity", "weather_pressure", "weather_precip",
  "hour", "is_weekend", "date", "functional_date",
  "Surgery_date", "Age", "Sex",
  "ema_missing"
)

df_ema <- df_ema_weather %>%
  select(all_of(ema_cols[ema_cols %in% names(df_ema_weather)])) %>%
  arrange(StudyID, time_block)

cat("\n  Final EMA lookup rows:", nrow(df_ema),
    "| Participants:", length(unique(df_ema$StudyID)), "\n")

# # -- 9d: missingness report - every column in df_ema ------------------------------
# cat("\n  df_ema column missingness:\n")
# ema_missingness <- tibble(
#   Column      = names(df_ema),
#   N_Missing   = sapply(df_ema, function(x) sum(is.na(x))),
#   N_Present   = sapply(df_ema, function(x) sum(!is.na(x))),
#   Missing_Pct = round(sapply(df_ema, function(x) mean(is.na(x)) * 100), 5)
# ) %>%
#   arrange(desc(Missing_Pct))
# 
# print(ema_missingness, n = nrow(ema_missingness))
# 
# write_csv(ema_missingness, file.path(report_dir, "df_ema_missingness.csv"))
# cat("  Missingness report saved to:",
#     file.path(report_dir, "df_ema_missingness.csv"), "\n")


# -- 9e: duplicate EMA filter ------------------------------------------------------
# When two completed surveys fall within 1 hour of each other for the same
# participant, the second (later) one is removed. The check is performed only
# on responded rows (!is.na(overall_pain)).
cat("\n  Filtering duplicate EMA responses within 1-hour windows...\n")

responded_rows <- df_ema %>%
  filter(!is.na(overall_pain)) %>%
  arrange(StudyID, time_block) %>%
  group_by(StudyID) %>%
  mutate(
    time_block_prior = lag(time_block, 1),
    time_since_prev  = as.numeric(
      difftime(time_block, time_block_prior, units = "mins")),
    is_duplicate = !is.na(time_since_prev) & time_since_prev < 60
  ) %>%
  ungroup()

# audit CSV - one row per removed duplicate
dup_audit <- responded_rows %>%
  filter(is_duplicate) %>%
  select(StudyID,
         time_block_removed = time_block,
         time_block_prior,
         gap_minutes        = time_since_prev) %>%
  arrange(StudyID, time_block_removed)

n_dups        <- nrow(dup_audit)
n_pids_dups   <- length(unique(dup_audit$StudyID))
dup_time_keys <- paste(dup_audit$StudyID, dup_audit$time_block_removed)

df_ema <- df_ema %>%
  mutate(.row_key = paste(StudyID, time_block)) %>%
  filter(!.row_key %in% dup_time_keys) %>%
  select(-.row_key)

cat(sprintf("  Removed %d duplicate responses across %d participants\n",
            n_dups, n_pids_dups))
cat("  Remaining EMA rows:", nrow(df_ema),
    "| Participants:", length(unique(df_ema$StudyID)), "\n")

write_csv(dup_audit, file.path(report_dir, "ema_duplicate_removal_audit.csv"))
cat("  Audit saved to:", file.path(report_dir, "ema_duplicate_removal_audit.csv"), "\n")

df_ema %>%
  mutate(time_block = format(time_block, "%Y-%m-%d %H:%M:%S")) %>%
  write_csv(ema_lookup_path)

cat("  Saved:", ema_lookup_path, "\n")


###############################################################################
# 10 - SUMMARY
###############################################################################

cat("\n=== Stage 1 complete ===\n")
cat("Reports written to:", report_dir, "\n")
cat("  df_expanded -", nrow(df_expanded), "rows,",
    length(unique(df_expanded$StudyID)), "participants\n")
cat("    ->", expanded_path, "\n")
cat("  df_ema      -", nrow(df_ema), "rows,",
    length(unique(df_ema$StudyID)), "participants\n")
cat("    ->", ema_lookup_path, "\n")
