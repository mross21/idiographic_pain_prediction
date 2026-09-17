# Idiographic Pain Prediction Analysis

This repository contains an R pipeline for building participant-specific models that predict whether pain will increase at the next EMA observation using Fitbit, sleep, EMA, and contextual data.

## Repository contents

- `install_packages.R`: installs the R packages used by the current codebase using the versions listed in the script
- `misc_tables.R`: produces supplemental participant demographics and predictor missingness tables for the model-eligible sample
- `mice.reuse.R`: local helper sourced by the modeling pipeline for reusing trained MICE imputations
- `sensitivity_analysis_data_window_minimum.R`: summarizes sensitivity analyses across valid-window minimum settings
- `step1_create_datasets.R`: builds the merged sensor and EMA modeling datasets
- `step2_imputation_and_modeling.R`: runs participant-specific imputation, feature generation, and modeling
- `step3_analysis.R`: aggregates participant-level results and produces summary outputs and figures

## Pipeline overview

The scripts are organized as a sequential workflow:

1. `step1_create_datasets.R`
   - Merges pre-op EMA composites with Fitbit heart rate/step data, sleep summaries, demographics, surgery dates, and daily weather.
   - Flags missing EMA, heart rate, steps, and sleep values.
   - Aggregates merged records into 5-minute blocks after removing likely Fitbit step artifacts.
   - Expands each participant to a complete 5-minute grid spanning their EMA study period.
   - Produces two key modeling inputs:
     - `df_expanded.csv`: full 5-minute sensor/context grid
     - `df_ema.csv`: EMA-level lookup table with current EMA values, context variables, and weather
   - Also writes QA reports on missingness and duplicate EMA removals.

2. `step2_imputation_and_modeling.R`
   - Loads `df_expanded.csv` and `df_ema.csv`.
   - Defines the prediction target as a binary indicator of whether `overall_pain` increased relative to the prior EMA.
   - Creates lagged EMA predictors, filters to participants with enough valid 1-hour pre-EMA heart-rate windows and enough sleep-covered days, and excludes participants with extreme class imbalance.
   - Uses rolling, forward-in-time folds within each participant.
   - Imputes missing Fitbit, sleep, and opioid data within each fold using `mice`/`mice.reuse` with mean-fill fallbacks.
   - Recomputes pre-EMA sensor features for each fold, then fits participant-specific models for four feature configurations:
     - Fitbit + context
     - Lagged EMA only
     - Fitbit + lagged EMA + context
     - Lag-1 pain only
   - Fits Random Forest, Elastic Net or Logistic Regression, Gaussian Process, and Ensemble models.
   - Saves participant-level predictions, pooled AUC and classification metrics, feature importance, fold caches, and summary tables.
   - Runs permuted null distributions for Lag-1 pain only feature configuration.

3. `step3_analysis.R`
   - Reads the participant-level outputs from step 2.
   - Aggregates model performance across participants and feature sets.
   - Summarizes top features from participant-best models.
   - Produces plots for pooled AUC distributions, participant-level performance, best-model frequencies, and feature importance.
   - Runs paired t-tests comparing pooled AUCs across analysis configurations.
   - Writes summary tables such as pooled performance metrics and sensitivity/specificity outputs.

## Modeling target

The current modeling script treats the outcome as a simple future pain escalation indicator:

- `Yes`: pain at the current EMA is greater than pain at the previous EMA
- `No`: pain does not increase

Prediction is evaluated idiographically, meaning each model is trained and tested separately within each participant over time.

## Setup

Before running the pipeline, install the required R packages with:

```r
source("install_packages.R")
```

The intended run order is:

1. `step1_create_datasets.R`
2. `step2_imputation_and_modeling.R`
3. `step3_analysis.R`

Supplemental analyses are provided in:

1. `misc_tables.R` (Demographics, predictor missingness)
2. `sensitivity_analysis_data_window_minimum.R` (Missingness by AUC, counts of AUC > 0.7/0.8 by feature configuration)
