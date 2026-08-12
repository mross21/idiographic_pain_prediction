################################################################################
#
# Personalized Pain Prediction - Step 2: Per-Participant Modeling
#
#   Loads the Stage 1 outputs (df_expanded, df_ema) written by
#   step1_create_datasets.R, determines which participants have enough sensor/EMA
#   coverage to model (window-based eligibility), builds per-fold imputed
#   feature sets, and trains per-participant models (Random Forest, Elastic
#   Net/Logistic Regression, Gaussian Process, and an ensemble) for each
#   analysis configuration.
#
# Run after step1_create_datasets.R. Requires df_expanded.csv and df_ema.csv.
################################################################################


###############################################################################
# 0 - LIBRARIES
###############################################################################

pkgs <- c("readr","dplyr","tidyr","lubridate","mice","caret",
          "randomForest","glmnet","kernlab","pROC","data.table","R.utils")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# mice.reuse is sourced locally because it may not be exported by miceadds.
source("mice.reuse.R")


###############################################################################
# 1 - PATHS + LOAD DATA
###############################################################################

output_dir     <- "/Users/f0085f6/Desktop/Frumkin_lab/personalizedPainPrediction_paper/PPP_Project/Output/pipeline/simple_pain_outcome"
fold_cache_dir <- file.path(output_dir, "folds")    # shared fold cache across threshold variants
results_dir    <- file.path(output_dir, "results")  # shared participant results across threshold variants

expanded_path   <- file.path(output_dir, "df_expanded.csv")
ema_lookup_path <- file.path(output_dir, "df_ema.csv")

if (!file.exists(expanded_path) || !file.exists(ema_lookup_path)) {
  stop("df_expanded.csv / df_ema.csv not found. Run step1_create_datasets.R first.", call. = FALSE)
}

cat("Loading Stage 1 outputs...\n")
df_expanded <- read_csv(expanded_path, show_col_types = FALSE) %>%
  mutate(time_block = as.POSIXct(time_block, tz = "UTC"))
df_ema <- read_csv(ema_lookup_path, show_col_types = FALSE) %>%
  mutate(time_block = as.POSIXct(time_block, tz = "UTC"))
cat("df_expanded:", nrow(df_expanded), "rows |",
    length(unique(df_expanded$StudyID)), "participants\n")
cat("df_ema:     ", nrow(df_ema),      "rows |",
    length(unique(df_ema$StudyID)),      "participants\n")

first_or_na <- function(x) {
  x <- na.omit(x)
  if (length(x) > 0) x[1] else NA_real_
}

write_log <- function(log_path, level, context, msg) {
  line <- sprintf("[%s] [%s] %s %s\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  level, context, msg)
  cat(line, file = log_path, append = TRUE)
}


###############################################################################
# 2 - ANALYSIS CONFIGURATIONS
###############################################################################

m_imp <- 5

# Eligibility is based on 1-hour pre-EMA sensor coverage rather than overall device missingness.
hr_window_min       <- 26                       # minimum valid HR windows
sleep_day_min       <- floor(hr_window_min / 5) # minimum EMA days with sleep data
analysis_label      <- paste0(hr_window_min, "ValidWindows")
analysis_output_dir <- file.path(output_dir, analysis_label)

fitbit_vars <- c(
  "averageHR","varianceHR",
  "totalSteps",
  "activeTime",
  "totalMinutesAsleep",
  "hr_missing","sleep_missing"
)

lag_ema_vars <- c(
  "overall_pain_lag1",  "overall_pain_lag2",
  "catastrophize_lag1", "catastrophize_lag2",
  "depress_lag1",       "depress_lag2",
  "interference_lag1",  "interference_lag2",
  "opioid_num_lag1"
)

context_vars <- c(
  "weather_temp","weather_humidity","weather_pressure","weather_precip",
  "hour","is_weekend", "days_since_first_ema"
)

analyses <- list(
  list(name  = "A1_FitbitOnly",
       vars  = c(fitbit_vars, context_vars),
       label = "Fitbit + Context"),
  list(name  = "A2_LagEMAOnly",
       vars  = lag_ema_vars,
       label = "Lag EMA Only"),
  list(name  = "A3_Combined",
       vars  = c(fitbit_vars, lag_ema_vars, context_vars),
       label = "Fitbit + Lag EMA + Context"),
  list(name  = "A4_Lag1PainOnly",
       vars  = "overall_pain_lag1",
       label = "Lag-1 Pain Only") 
)

if (!dir.exists(analysis_output_dir))
  dir.create(analysis_output_dir, recursive = TRUE)
if (!dir.exists(fold_cache_dir))
  dir.create(fold_cache_dir, recursive = TRUE)
if (!dir.exists(results_dir))
  dir.create(results_dir, recursive = TRUE)

# Logs are shared across threshold variants.
log_imputation <- file.path(fold_cache_dir, "log_imputation.txt")
log_modeling   <- file.path(results_dir,    "log_modeling.txt")

for (a in analyses) {
  d_results <- file.path(results_dir, a$name, "PatientResults")
  if (!dir.exists(d_results)) dir.create(d_results, recursive = TRUE)
  d_summary <- file.path(analysis_output_dir, a$name, "Summary")
  if (!dir.exists(d_summary)) dir.create(d_summary, recursive = TRUE)
}

###############################################################################
# 3 - HELPER FUNCTIONS
###############################################################################

compute_ema_features <- function(pid, t, sensor_df, window_hours = 1) {
  na_row <- function() data.frame(
    StudyID = pid, time_block = t,
    averageHR = NA_real_, varianceHR = NA_real_, wearTime = NA_real_,
    totalSteps = NA_real_, averageSteps = NA_real_, varianceSteps = NA_real_,
    sedentaryTime = NA_real_, activeTime = NA_real_,
    totalMinutesAsleep = NA_real_,
    hr_missing = NA_real_, steps_missing = NA_real_, sleep_missing = NA_real_
  )
  tMin <- t - hours(window_hours)
  sGrp <- sensor_df %>% filter(time_block > tMin & time_block <= t)
  if (sum(!is.na(sGrp$hr) & !is.na(sGrp$steps)) < 2) return(na_row())
  avgHR <- mean(sGrp$hr, na.rm = TRUE)
  varHR <- var(sGrp$hr,  na.rm = TRUE)
  dt           <- diff(sGrp$time_block)
  n_consec     <- sum(as.numeric(dt, units = "mins") == 5)
  wearTime     <- n_consec * 5
  if (wearTime == 0) return(na_row())
  totalSteps      <- sum(sGrp$total_steps,  na.rm = TRUE)
  avgSteps        <- mean(sGrp$steps, na.rm = TRUE)
  varSteps        <- var(sGrp$steps,  na.rm = TRUE)
  fbSed  <- sGrp %>% filter(steps == 0)
  sedTime <- sum(as.numeric(diff(fbSed$time_block), units = "mins") == 5) * 5
  fbAct  <- sGrp %>% filter(steps > 50)
  actTime <- sum(as.numeric(diff(fbAct$time_block), units = "mins") == 5) * 5
  sleep_val <- sGrp %>% filter(time_block == max(time_block)) %>%
    pull(totalMinutesAsleep) %>% first_or_na()
  data.frame(
    StudyID = pid, time_block = t,
    averageHR = avgHR, varianceHR = varHR, wearTime = wearTime,
    totalSteps = totalSteps, averageSteps = avgSteps, varianceSteps = varSteps,
    sedentaryTime = sedTime, activeTime = actTime,
    totalMinutesAsleep = sleep_val,
    hr_missing    = mean(sGrp$hr_missing,    na.rm = TRUE),
    steps_missing = mean(sGrp$steps_missing, na.rm = TRUE),
    sleep_missing = round(mean(sGrp$sleep_missing, na.rm = TRUE)) # get majority vote
  )
}

pool_auc_rubins <- function(aucs, ses) {
  m    <- length(aucs)
  eps  <- 1e-6
  aucs <- pmax(eps, pmin(1 - eps, aucs))
  Z    <- log(aucs / (1 - aucs))
  SE_Z <- ses / (aucs * (1 - aucs))
  Z_bar <- mean(Z)
  W     <- mean(SE_Z^2)
  B     <- var(Z)
  T_var <- W + (1 + 1/m) * B
  SE    <- sqrt(T_var)
  inv_f <- function(z) 1 / (1 + exp(-z))
  list(auc = inv_f(Z_bar),
       ci_lo = inv_f(Z_bar - 1.96 * SE),
       ci_hi = inv_f(Z_bar + 1.96 * SE),
       se = SE, within_var = W, between_var = B)
}

pool_proportion_rubins <- function(vals, ns = NULL) {
  m    <- length(vals)
  eps  <- 1e-6
  vals <- pmax(eps, pmin(1 - eps, vals))
  L    <- log(vals / (1 - vals))
  L_bar <- mean(L)
  W     <- if (!is.null(ns)) mean(1 / (ns * vals * (1 - vals))) else
    mean(1 / (50 * vals * (1 - vals)))
  B     <- var(L)
  T_var <- W + (1 + 1/m) * B
  SE    <- sqrt(T_var)
  inv   <- function(x) 1 / (1 + exp(-x))
  list(est   = inv(L_bar),
       ci_lo = inv(L_bar - 1.96 * SE),
       ci_hi = inv(L_bar + 1.96 * SE))
}

pool_kappa_rubins <- function(kappas, kappa_ses) {
  m      <- length(kappas)
  eps    <- 1e-6
  kappas <- pmax(-(1 - eps), pmin(1 - eps, kappas))
  Z      <- 0.5 * log((1 + kappas) / (1 - kappas))
  SE_Z   <- kappa_ses / (1 - kappas^2)
  Z_bar  <- mean(Z)
  W      <- mean(SE_Z^2)
  B      <- var(Z)
  T_var  <- W + (1 + 1/m) * B
  SE     <- sqrt(T_var)
  inv_k  <- function(z) (exp(2*z) - 1) / (exp(2*z) + 1)
  list(est      = inv_k(Z_bar),
       ci_lo    = inv_k(Z_bar - 1.96 * SE),
       ci_hi    = inv_k(Z_bar + 1.96 * SE),
       se       = SE,
       within_var  = W,
       between_var = B)
}

create_rolling_folds <- function(data, min_train = 20, val_size = 10,
                                 test_size = 1, expanding = TRUE) {
  n <- nrow(data)
  folds <- list()
  fc <- 0
  for (ts in (min_train + val_size + 1):(n - test_size + 1)) {
    vs <- ts - val_size
    tr <- if (expanding) 1:(vs - 1) else max(1, vs - min_train):(vs - 1)
    if (length(tr) >= min_train) {
      fc <- fc + 1
      folds[[fc]] <- list(train      = tr,
                          validation = vs:(ts - 1),
                          test       = ts:(ts + test_size - 1),
                          fold_id    = fc)
    }
  }
  folds
}

convert_predictor_types <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    if (is.logical(df[[v]]))   df[[v]] <- as.numeric(df[[v]])
    if (is.character(df[[v]])) df[[v]] <- as.numeric(as.factor(df[[v]]))
  }
  df
}

clean_participant_data <- function(pdata, pvars) {
  valid <- pvars[pvars %in% names(pdata)]
  for (v in valid) {
    if (is.numeric(pdata[[v]])) {
      uv <- unique(pdata[[v]][!is.na(pdata[[v]])])
      if (length(uv) <= 1 ||
          max(table(pdata[[v]], useNA = "no")) / sum(!is.na(pdata[[v]])) > 0.95)
        valid <- setdiff(valid, v)
    }
  }
  list(data = pdata, valid_predictors = valid)
}

calculate_permutation_importance <- function(model, X_test, y_test,
                                             n_perm = 10, type = "rf") {
  set.seed(126)
  X_test <- as.data.frame(X_test)
  get_p <- function(m, X) tryCatch({
    if (type == "rf")          predict(m, X, type = "prob")$Yes
    else if (type == "glmnet") predict(m, as.matrix(X), type = "prob")$Yes
    else if (type == "glm")    predict(m, X, type = "prob")$Yes
    else if (type == "gp")     predict(m, X, type = "prob")$Yes
    else rep(0.5, nrow(X))
  }, error = function(e) rep(0.5, nrow(X)))
  safe_auc <- function(resp, pred) tryCatch(
    as.numeric(auc(roc(resp, pred, quiet = TRUE))), error = function(e) 0.5)
  base   <- safe_auc(y_test, get_p(model, X_test))
  scores <- setNames(numeric(ncol(X_test)), names(X_test))
  for (feat in names(scores)) {
    pa <- numeric(n_perm)
    for (pp in seq_len(n_perm)) {
      set.seed(126 + pp); Xp <- X_test; Xp[[feat]] <- sample(Xp[[feat]])
      pa[pp] <- safe_auc(y_test, get_p(model, Xp))
    }
    scores[feat] <- max(0, base - mean(pa))
  }
  scores * 100
}


###############################################################################
# 4 - PER-PARTICIPANT MODEL BUILDER
###############################################################################

build_participant_model <- function(pid, predictor_vars, patient_dir,
                                    analysis_name,
                                    fold_rds_files,
                                    holdout_rds_file = NULL) {
  if (length(fold_rds_files) == 0) {
    message("  Skipping ", pid, ": no fold RDS files")
    return(NULL)
  }
  
  last_fold        <- readRDS(fold_rds_files[[length(fold_rds_files)]])
  last_train       <- convert_predictor_types(last_fold[[1]]$train, predictor_vars)
  cleaned          <- clean_participant_data(last_train, predictor_vars)
  valid_predictors <- cleaned$valid_predictors
  rm(last_fold, last_train)
  
  if (length(valid_predictors) < 1) {
    message("  Skipping ", pid, ": no valid predictors")
    return(NULL)
  }
  
  phash     <- sum(as.numeric(charToRaw(pid)))
  base_seed <- 126 + phash
  
  all_imp_preds <- vector("list", m_imp)
  stored_models <- vector("list", m_imp)
  for (ii in seq_len(m_imp)) {
    all_imp_preds[[ii]] <- data.frame()
    stored_models[[ii]] <- list()
  }
  
  is_a4 <- (analysis_name == "A4_Lag1PainOnly")
  
  for (fi in seq_along(fold_rds_files)) {
    fold_data <- readRDS(fold_rds_files[[fi]])
    
    for (imp_idx in seq_len(m_imp)) {
      if (imp_idx > length(fold_data)) next
      imp_split <- fold_data[[imp_idx]]
      train_d   <- convert_predictor_types(imp_split$train, valid_predictors)
      val_d     <- convert_predictor_types(imp_split$val,   valid_predictors)
      test_d    <- convert_predictor_types(imp_split$test,  valid_predictors)
      
      base_prob <- mean(train_d$pain_increasing == "Yes", na.rm = TRUE)
      if (min(table(train_d$pain_increasing)) < 2) next
      
      y_train <- train_d$pain_increasing
      X_train <- as.data.frame(train_d %>% select(all_of(valid_predictors)))
      y_val   <- val_d$pain_increasing
      X_val   <- as.data.frame(val_d   %>% select(all_of(valid_predictors)))
      y_test  <- test_d$pain_increasing
      X_test  <- as.data.frame(test_d  %>% select(all_of(valid_predictors)))
      
      const_v <- sapply(X_train, function(x) length(unique(x)) <= 1)
      if (sum(!const_v) < 1) next
      if (any(const_v)) {
        X_train <- X_train[, !const_v, drop = FALSE]
        X_val   <- X_val[,   !const_v, drop = FALSE]
        X_test  <- X_test[,  !const_v, drop = FALSE]
      }
      fold_predictors <- names(X_train)
      
      set.seed(base_seed + (imp_idx - 1) * 1000 + fi)
      
      # Inner tuning uses forward-in-time resampling within each participant.
      chunk_size <- max(5L, floor(nrow(X_train) / 4))
      inner_cv   <- trainControl(
        method          = "timeslice",
        initialWindow   = chunk_size,
        horizon         = chunk_size,
        fixedWindow     = FALSE,
        skip            = 1,
        classProbs      = TRUE,
        summaryFunction = twoClassSummary,
        allowParallel   = FALSE
      )
      
      ctx_fold <- sprintf("[%s] [pid=%s] [fold=%02d] [imp=%d]",
                          analysis_name, pid, fi, imp_idx)
      fold_result <- withCallingHandlers(
        tryCatch({
          # Random Forest
          rf_mod <- tryCatch(
            R.utils::withTimeout(
              train(x = X_train, y = y_train, method = "rf",
                    tuneGrid  = data.frame(mtry = unique(c(
                      max(1L, floor(sqrt(ncol(X_train)))),
                      max(1L, floor(ncol(X_train) / 2))
                    ))),
                    metric    = "ROC", ntree = 200,
                    trControl = inner_cv),
              timeout = 30, # seconds
              onTimeout = "silent"
            ),
            error = function(e) {
              message("RF model failed or timed out: ", conditionMessage(e))
              NULL
            }
          )
          rf_val  <- predict(rf_mod, X_val,  type = "prob")$Yes
          rf_test <- predict(rf_mod, X_test, type = "prob")$Yes
          
          sc  <- scale(X_train)
          ctr <- attr(sc, "scaled:center")
          scl <- attr(sc, "scaled:scale"); scl[scl == 0] <- 1
          Xvs <- as.data.frame(scale(X_val,  center = ctr, scale = scl))
          Xts <- as.data.frame(scale(X_test, center = ctr, scale = scl))
          sc  <- as.data.frame(sc)
          
          # Elastic Net
          en_mod <- if (is_a4) {
            tryCatch(
              R.utils::withTimeout(
                train(x = sc, y = y_train, method = "glm",
                      metric    = "ROC",
                      trControl = trainControl(method = "none", classProbs = TRUE,
                                               summaryFunction = twoClassSummary),
                      family    = "binomial"),
                timeout = 30, # seconds
                onTimeout = "silent"
              ),
              error = function(e) {
                message("LR model failed or timed out: ", conditionMessage(e))
                NULL
              }
            )
          } else if (all(table(y_train) >= 8)) {
            tryCatch(
              R.utils::withTimeout(
                train(x = as.matrix(sc), y = y_train, method = "glmnet",
                      tuneGrid  = expand.grid(alpha  = c(0.5, 1.0),
                                              lambda = c(0.01, 0.1, 1.0)),
                      metric    = "ROC",
                      trControl = inner_cv),
                timeout = 30, # seconds
                onTimeout = "silent"
              ),
              error = function(e) {
                message("EN model failed or timed out: ", conditionMessage(e))
                NULL
              }
            )
          } else NULL
          
          predict_en <- function(Xs) {
            if (is.null(en_mod)) return(rep(base_prob, nrow(Xs)))
            if (is_a4) predict(en_mod, Xs,             type = "prob")$Yes
            else       predict(en_mod, as.matrix(Xs), type = "prob")$Yes
          }
          en_val  <- predict_en(Xvs)
          en_test <- predict_en(Xts)
          
          # Gaussian Process
          sigma_est <- kernlab::sigest(as.matrix(sc), scaled = FALSE)
          gp_mod <- tryCatch(
            R.utils::withTimeout(
              train(x = sc, y = y_train, method = "gaussprRadial",
                    tuneGrid  = data.frame(sigma = c(
                      sigma_est[1], sigma_est[2], sigma_est[3])),
                    metric    = "ROC",
                    trControl = inner_cv),
              timeout = 30, # seconds
              onTimeout = "silent"
            ),
            error = function(e) {
              message("GP model failed or timed out: ", conditionMessage(e))
              NULL
            }
          )
          
          gp_val  <- if (!is.null(gp_mod))
            predict(gp_mod, Xvs, type = "prob")$Yes else
              rep(base_prob, nrow(X_val))
          gp_test <- if (!is.null(gp_mod))
            predict(gp_mod, Xts, type = "prob")$Yes else
              rep(base_prob, nrow(X_test))
          
          # Ensemble
          meta_y <- as.numeric(y_val == "Yes")
          meta_X <- data.frame(RF = rf_val, EN = en_val, GP = gp_val)
          ens_test <- tryCatch({
            if (!any(is.na(meta_X)) && length(unique(meta_y)) >= 2 &&
                min(table(meta_y)) >= 3) {
              mm <- suppressWarnings(
                glm(meta_y ~ RF + EN + GP, data = meta_X,
                    family  = binomial(),
                    control = glm.control(maxit = 100, epsilon = 1e-6))
              )
              fitted_p <- fitted(mm)
              if (any(fitted_p <= 0.001 | fitted_p >= 0.999)) {
                0.4*rf_test + 0.4*en_test + 0.2*gp_test
              } else {
                predict(mm, data.frame(RF = rf_test, EN = en_test, GP = gp_test),
                        type = "response")
              }
            } else 0.4*rf_test + 0.4*en_test + 0.2*gp_test
          }, error = function(e) 0.4*rf_test + 0.4*en_test + 0.2*gp_test)
          
          en_pred_col <- if (is_a4) "LR_Pred" else "EN_Pred"
          fold_res <- data.frame(
            Date       = test_d$functional_date,
            Observed   = as.integer(y_test == "Yes"),
            RF_Pred    = as.numeric(rf_test),
            EN_or_LR   = as.numeric(en_test),
            GP_Pred    = as.numeric(gp_test),
            Ens_Pred   = as.numeric(ens_test),
            BaseProb   = base_prob,
            Fold       = fi,
            TrainSize  = nrow(train_d),
            Imputation = imp_idx
          )
          names(fold_res)[names(fold_res) == "EN_or_LR"] <- en_pred_col
          list(
            fold_res   = fold_res,
            fold_model = list(
              rf = rf_mod, glmnet = en_mod, gp = gp_mod,
              scaling = list(center = ctr, scale = scl),
              predictors = fold_predictors,
              is_logistic = is_a4)
          )
        }, error = function(e) {
          message("  Error fold ", fi, " imp ", imp_idx, ": ", e$message)
          write_log(log_modeling, "ERROR", ctx_fold, e$message)
          NULL
        }),
        warning = function(w) {
          write_log(log_modeling, "WARN", ctx_fold, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          write_log(log_modeling, "MSG", ctx_fold, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      if (!is.null(fold_result)) {
        all_imp_preds[[imp_idx]] <- bind_rows(all_imp_preds[[imp_idx]], fold_result$fold_res)
        stored_models[[imp_idx]][[sprintf("fold_%04d", fi)]] <- fold_result$fold_model
      }
    }
  }
  
  valid_imps <- which(
    sapply(all_imp_preds, function(x) is.data.frame(x) && nrow(x) > 0))
  if (length(valid_imps) == 0) {
    message("  No valid imputation results for ", pid); return(NULL)
  }
  
  
  model_cols        <- c("RF_Pred", if (is_a4) "LR_Pred" else "EN_Pred", "GP_Pred","Ens_Pred")
  model_names_clean <- if (is_a4)
    c("Random Forest","Logistic Regression","Gaussian Process","Ensemble") else
      c("Random Forest","Elastic Net","Gaussian Process","Ensemble")
  
  pooled_aucs <- setNames(vector("list", 4), model_names_clean)
  for (mc in seq_along(model_cols)) {
    aucs <- numeric(length(valid_imps))
    ses  <- numeric(length(valid_imps))
    for (ii in seq_along(valid_imps)) {
      pd <- all_imp_preds[[valid_imps[ii]]]
      tryCatch({
        r         <- roc(pd$Observed, pd[[model_cols[mc]]], quiet = TRUE)
        aucs[ii]  <- as.numeric(auc(r))
        ses[ii]   <- sqrt(var(r))
      }, error = function(e) { aucs[ii] <<- 0.5; ses[ii] <<- 0.05 })
    }
    pooled_aucs[[model_names_clean[mc]]] <- pool_auc_rubins(aucs, ses)
  }
  
  calc_metrics_one_imp <- function(obs, pred) {
    tryCatch({
      if (length(unique(obs)) < 2 || length(unique(pred)) < 2)
        return(list(auc=0.5, sens=NA, spec=NA, ppv=NA, npv=NA,
                    f1=NA, kappa=NA, kappa_se=NA, n=length(obs),
                    npos=sum(obs==1), nneg=sum(obs==0)))
      r  <- roc(obs, pred, quiet = TRUE)
      co <- coords(r, "best", ret = c("threshold","sensitivity","specificity"),
                   best.method = "youden")
      co <- co[1, , drop = FALSE]
      pc <- ifelse(pred >= co$threshold, 1, 0)
      TP <- sum(pc==1&obs==1); TN <- sum(pc==0&obs==0)
      FP <- sum(pc==1&obs==0); FN <- sum(pc==0&obs==1)
      ppv <- ifelse(TP+FP>0, TP/(TP+FP), NA)
      npv <- ifelse(TN+FN>0, TN/(TN+FN), NA)
      f1  <- ifelse(!is.na(ppv) && co$sensitivity+ppv>0,
                    2*ppv*co$sensitivity/(ppv+co$sensitivity), NA)
      n   <- length(obs)
      py  <- ((TP+FN)/n)*((TP+FP)/n); pn <- ((TN+FP)/n)*((TN+FN)/n)
      kp  <- ifelse(py+pn<1, ((TP+TN)/n-(py+pn))/(1-(py+pn)), NA)
      pe      <- py + pn
      kp_se   <- ifelse(!is.na(kp) && pe < 1,
                        sqrt(pe * (1 - pe) / (n * (1 - pe)^2)),
                        NA_real_)
      list(auc=as.numeric(auc(r)), sens=co$sensitivity, spec=co$specificity,
           ppv=ppv, npv=npv, f1=f1, kappa=kp, kappa_se=kp_se,
           n=n, npos=sum(obs==1), nneg=sum(obs==0))
    }, error = function(e)
      list(auc=0.5, sens=NA, spec=NA, ppv=NA, npv=NA,
           f1=NA, kappa=NA, kappa_se=NA, n=length(obs), npos=NA, nneg=NA))
  }
  
  pool_one_model <- function(col, model_name) {
    mets <- lapply(valid_imps, function(ii)
      calc_metrics_one_imp(all_imp_preds[[ii]]$Observed, all_imp_preds[[ii]][[col]]))
    ns_all <- sapply(mets, function(x) x$n)
    pm <- function(key) {
      v  <- sapply(mets, function(x) x[[key]])
      ok <- !is.na(v)
      v  <- v[ok]
      if (length(v) == 0) return(list(est=NA, ci_lo=NA, ci_hi=NA))
      pool_proportion_rubins(v, ns_all[ok])
    }
    auc_p <- pool_auc_rubins(sapply(mets, function(x) x$auc),
                             rep(0.05, length(mets)))
    kappa_p <- {
      kv  <- sapply(mets, function(x) x$kappa)
      ksv <- sapply(mets, function(x) x$kappa_se)
      ok  <- !is.na(kv) & !is.na(ksv)
      if (sum(ok) == 0) list(est=NA_real_, ci_lo=NA_real_, ci_hi=NA_real_)
      else pool_kappa_rubins(kv[ok], ksv[ok])
    }
    data.frame(
      StudyID     = pid, Model = model_name,
      AUC         = auc_p$auc,  AUC_CI_lo  = auc_p$ci_lo,  AUC_CI_hi  = auc_p$ci_hi,
      Sensitivity = pm("sens")$est, Sens_CI_lo = pm("sens")$ci_lo, Sens_CI_hi = pm("sens")$ci_hi,
      Specificity = pm("spec")$est, Spec_CI_lo = pm("spec")$ci_lo, Spec_CI_hi = pm("spec")$ci_hi,
      PPV         = pm("ppv")$est,  PPV_CI_lo  = pm("ppv")$ci_lo,  PPV_CI_hi  = pm("ppv")$ci_hi,
      NPV         = pm("npv")$est,  NPV_CI_lo  = pm("npv")$ci_lo,  NPV_CI_hi  = pm("npv")$ci_hi,
      F1          = pm("f1")$est,   F1_CI_lo   = pm("f1")$ci_lo,   F1_CI_hi   = pm("f1")$ci_hi,
      Kappa       = kappa_p$est,  Kappa_CI_lo = kappa_p$ci_lo,  Kappa_CI_hi = kappa_p$ci_hi,
      N_Obs = mets[[1]]$n, N_Pos = mets[[1]]$npos, N_Neg = mets[[1]]$nneg
    )
  }
  participant_sens_spec <- bind_rows(
    mapply(pool_one_model, model_cols, model_names_clean, SIMPLIFY = FALSE))
  
  aggregated_importance <- data.frame()
  en_label <- if (is_a4) "Logistic Regression" else "Elastic Net"
  
  imp_scores_by_model <- setNames(
    list(list(), list(), list()),
    c("Random Forest", en_label, "Gaussian Process")
  )
  
  # Hold out the final five complete EMA observations for permutation importance only.
  ema_holdout <- if (!is.null(holdout_rds_file) && file.exists(holdout_rds_file))
    readRDS(holdout_rds_file) else NULL
  
  if (is.null(ema_holdout) || length(ema_holdout) == 0) {
    cat("  Skipping permutation importance - holdout set unavailable\n")
  } else {
    for (ii in valid_imps) {
      fold_models_ii <- stored_models[[ii]]
      last_fold_key  <- tail(names(fold_models_ii), 1)
      last_fold_ii   <- fold_models_ii[[last_fold_key]]
      if (is.null(last_fold_ii)) next
      
      # use the imputation-matched holdout split
      holdout_ii <- if (ii <= length(ema_holdout)) ema_holdout[[ii]] else ema_holdout[[1]]
      if (is.null(holdout_ii) || nrow(holdout_ii) < 2) next
      if (length(unique(holdout_ii$pain_increasing)) < 2) next
      
      # align X_ii columns to exactly match the model's predictor set
      imp_preds_ii <- if (!is.null(last_fold_ii$predictors))
        last_fold_ii$predictors else valid_predictors
      X_ii <- holdout_ii %>%
        select(any_of(imp_preds_ii)) %>%
        as.data.frame()
      for (mc in setdiff(imp_preds_ii, names(X_ii))) X_ii[[mc]] <- NA_real_
      X_ii <- X_ii[, imp_preds_ii, drop = FALSE]
      X_ii <- X_ii %>% filter(rowSums(is.na(.)) == 0)
      y_ii <- as.integer(holdout_ii$pain_increasing[
        rowSums(is.na(holdout_ii %>% select(any_of(imp_preds_ii)))) == 0] == "Yes")
      
      if (nrow(X_ii) < 2 || length(unique(y_ii)) < 2) next
      
      if (!is.null(last_fold_ii$rf)) {
        ri <- tryCatch(
          calculate_permutation_importance(last_fold_ii$rf, X_ii, y_ii, type="rf"),
          error = function(e) { cat("  [ERROR RF importance]:", e$message, "\n"); NULL }
        )
        if (!is.null(ri))
          imp_scores_by_model[["Random Forest"]] <-
            c(imp_scores_by_model[["Random Forest"]], list(ri))
      }
      if (!is.null(last_fold_ii$glmnet)) {
        if (isTRUE(last_fold_ii$is_logistic)) {
          ei <- tryCatch(
            calculate_permutation_importance(last_fold_ii$glmnet, X_ii, y_ii, type="glm"),
            error = function(e) { cat("  [ERROR EN importance]:", e$message, "\n"); NULL }
          )
        } else {
          scale_vars <- names(last_fold_ii$scaling$center)
          Xs_ii <- X_ii[, intersect(scale_vars, names(X_ii)), drop = FALSE]
          for (mc in setdiff(scale_vars, names(X_ii))) Xs_ii[[mc]] <- NA_real_
          Xs_ii <- as.data.frame(scale(Xs_ii[, scale_vars, drop=FALSE],
                                       center = last_fold_ii$scaling$center,
                                       scale  = last_fold_ii$scaling$scale))
          ei <- tryCatch(
            calculate_permutation_importance(last_fold_ii$glmnet, Xs_ii, y_ii, type="glmnet"),
            error = function(e) { cat("  [ERROR EN importance]:", e$message, "\n"); NULL }
          )
        }
        if (!is.null(ei))
          imp_scores_by_model[[en_label]] <-
            c(imp_scores_by_model[[en_label]], list(ei))
      }
      if (!is.null(last_fold_ii$gp)) {
        scale_vars <- names(last_fold_ii$scaling$center)
        Xs_ii <- X_ii[, intersect(scale_vars, names(X_ii)), drop = FALSE]
        for (mc in setdiff(scale_vars, names(X_ii))) Xs_ii[[mc]] <- NA_real_
        Xs_ii <- as.data.frame(scale(Xs_ii[, scale_vars, drop=FALSE],
                                     center = last_fold_ii$scaling$center,
                                     scale  = last_fold_ii$scaling$scale))
        gi <- tryCatch(
          calculate_permutation_importance(last_fold_ii$gp, Xs_ii, y_ii, type="gp"),
          error = function(e) { cat("  [ERROR GP importance]:", e$message, "\n"); NULL }
        )
        if (!is.null(gi))
          imp_scores_by_model[["Gaussian Process"]] <-
          c(imp_scores_by_model[["Gaussian Process"]], list(gi))
      }
    }
  }
  
  imp_parts <- list()
  for (mn in names(imp_scores_by_model)) {
    scores_list <- imp_scores_by_model[[mn]]
    if (length(scores_list) == 0) next
    all_vars   <- names(scores_list[[1]])
    avg_scores <- rowMeans(
      do.call(cbind, lapply(scores_list, function(s) s[all_vars])),
      na.rm = TRUE)
    imp_parts[[mn]] <- data.frame(
      Model = mn, Variable = all_vars, Mean_Importance = as.numeric(avg_scores))
  }
  
  if (length(imp_parts) > 0) {
    aggregated_importance <- bind_rows(imp_parts)
    ens_imp <- aggregated_importance %>%
      group_by(Variable) %>%
      summarise(Mean_Importance = mean(Mean_Importance, na.rm=TRUE), .groups="drop") %>%
      mutate(Model = "Ensemble")
    aggregated_importance <- bind_rows(aggregated_importance, ens_imp)
  }
  
  pdir <- file.path(patient_dir, pid)
  if (!dir.exists(pdir)) dir.create(pdir)
  
  bind_rows(all_imp_preds[valid_imps]) %>%
    write_csv(file.path(pdir, paste0(pid,"_predictions_all_imputations.csv")))
  
  bind_rows(lapply(seq_along(valid_imps), function(ii) {
    pd  <- all_imp_preds[[valid_imps[ii]]]
    row <- data.frame(Imputation = valid_imps[ii])
    for (mc in seq_along(model_cols)) tryCatch({
      row[[model_names_clean[mc]]] <-
        as.numeric(auc(roc(pd$Observed, pd[[model_cols[mc]]], quiet=TRUE)))
    }, error = function(e) { row[[model_names_clean[mc]]] <<- 0.5 })
    row
  })) %>%
    write_csv(file.path(pdir, paste0(pid,"_auc_per_imputation.csv")))
  
  bind_rows(lapply(model_names_clean, function(mn) {
    pa <- pooled_aucs[[mn]]
    data.frame(Model=mn, Pooled_AUC=pa$auc,
               CI_lo=pa$ci_lo, CI_hi=pa$ci_hi,
               Within_Var=pa$within_var, Between_Var=pa$between_var)
  })) %>%
    write_csv(file.path(pdir, paste0(pid,"_pooled_auc.csv")))
  
  participant_sens_spec %>%
    write_csv(file.path(pdir, paste0(pid,"_pooled_sens_spec.csv")))
  
  if (nrow(aggregated_importance) > 0)
    write_csv(aggregated_importance,
              file.path(pdir, paste0(pid,"_feature_importance.csv")))
  
  en_name   <- model_names_clean[2]
  en_col    <- if (is_a4) "LogisticRegression_AUC" else "ElasticNet_AUC"
  en_ci_col <- if (is_a4) "LR_CI" else "EN_CI"
  best_base <- max(pooled_aucs[["Random Forest"]]$auc,
                   pooled_aucs[[en_name]]$auc,
                   pooled_aucs[["Gaussian Process"]]$auc)
  result <- data.frame(
    StudyID             = pid,
    RandomForest_AUC    = pooled_aucs[["Random Forest"]]$auc,
    RF_CI               = paste0("[",round(pooled_aucs[["Random Forest"]]$ci_lo,3),
                                 ", ",round(pooled_aucs[["Random Forest"]]$ci_hi,3),"]"),
    stringsAsFactors    = FALSE
  )
  result[[en_col]]    <- pooled_aucs[[en_name]]$auc
  result[[en_ci_col]] <- paste0("[",round(pooled_aucs[[en_name]]$ci_lo,3),
                                ", ",round(pooled_aucs[[en_name]]$ci_hi,3),"]")
  result$GaussianProcess_AUC = pooled_aucs[["Gaussian Process"]]$auc
  result$GP_CI               = paste0("[",round(pooled_aucs[["Gaussian Process"]]$ci_lo,3),
                                      ", ",round(pooled_aucs[["Gaussian Process"]]$ci_hi,3),"]")
  result$Ensemble_AUC        = pooled_aucs[["Ensemble"]]$auc
  result$Ens_CI              = paste0("[",round(pooled_aucs[["Ensemble"]]$ci_lo,3),
                                      ", ",round(pooled_aucs[["Ensemble"]]$ci_hi,3),"]")
  result$BaselineAUC         = 0.5
  result$Best_Base_AUC       = best_base
  result$Ensemble_Improvement = pooled_aucs[["Ensemble"]]$auc - best_base
  result$NumObservations     = nrow(all_imp_preds[[valid_imps[1]]])
  result$N_Imputations       = length(valid_imps)
  
  saveRDS(result, file.path(pdir, paste0(pid, "_participant_summary_row.rds")))
  
  result
}


###############################################################################
# 5 - PRE-LOOP: ELIGIBILITY (5a) + FOLD-LEVEL IMPUTATION + FEATURES (5b)
###############################################################################

cat("\n=== Computing eligibility (observed sequence) ===\n")

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
    # pain_flag         = case_when(
    #   is.na(overall_pain_lag1) | is.na(diff) ~ NA_integer_,
    #   (overall_pain_lag1 >= 30 & diff >= 0)  |
    #     (overall_pain_lag1 <  30 & diff >= 10) |
    #     (overall_pain_lag1 == 0  & diff >= 10) ~ 1L,
    #   TRUE ~ 0L),
    total_transition  = max(cumsum(replace(pain_flag, is.na(pain_flag), 0))),
    n                 = n(),
    perc_transition   = total_transition / n,
    pain_increasing   = factor(
      case_when(
        is.na(pain_flag) ~ NA_character_,
        pain_flag == 1L  ~ "Yes",
        TRUE             ~ "No"
      ), levels = c("No","Yes")),
    # Compute EMA lags before filtering so gaps do not reset the sequence.
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

cat("Class distribution (observed sequence):",
    table(base_data$pain_increasing), "\n")

lag_cols_all <- c("overall_pain_lag1",
                  "overall_pain_lag2",
                  "catastrophize_lag1", "catastrophize_lag2",
                  "depress_lag1",       "depress_lag2",
                  "interference_lag1",  "interference_lag2",
                  "ema_missing_lag1")

base_responded <- base_data %>% filter(!is.na(pain_increasing))
total_rows     <- nrow(base_responded)

na_report <- bind_rows(lapply(lag_cols_all, function(col) {
  n_na <- sum(is.na(base_responded[[col]]))
  data.frame(
    EMA_Item     = col,
    NA_rows      = n_na,
    Total_rows   = total_rows,
    Pct_removed  = round(n_na / total_rows * 100, 2)
  )
}))

# also report pain_increasing NAs - base_data is already filtered to !is.na(overall_pain)
n_na_pi <- sum(is.na(base_data$pain_increasing))
na_report <- bind_rows(
  data.frame(EMA_Item    = "pain_increasing",
             NA_rows     = n_na_pi,
             Total_rows  = nrow(base_data),
             Pct_removed = round(n_na_pi / nrow(base_data) * 100, 2)),
  na_report
)

write_csv(na_report, file.path(output_dir, "ema_lag_na_removal_report.csv"))
cat("EMA lag NA removal report:\n")
print(na_report)

window_missingness_csv <- file.path(output_dir, "window_missingness_summary.csv")

if (file.exists(window_missingness_csv)) {
  cat("  Loading cached window missingness summary...\n")
  window_eligibility <- read_csv(window_missingness_csv, show_col_types = FALSE)
} else {
  cat("\n  Computing window-based sensor missingness for eligibility...\n")
  
  ema_timestamps <- base_data %>%
    filter(!is.na(pain_increasing),
           !is.na(catastrophize),
           !is.na(depress),
           !is.na(interference)) %>%
    select(StudyID, time_block, functional_date) %>%
    mutate(window_start = time_block - hours(1),
           window_end   = time_block)
  
  # vectorised non-equi join: assign each sensor row to its EMA window in one pass
  # this replaces a row-by-row filter loop (O(n_ema * n_sensor) -> O(n_sensor log n))
  sensor_dt  <- as.data.table(df_expanded)[,
                                           .(StudyID, time_block, hr_missing, sleep_missing)]
  windows_dt <- as.data.table(ema_timestamps)
  
  # non-equi join: sensor row falls in window if time_block > window_start & <= window_end
  matched <- sensor_dt[windows_dt,
                       on = .(StudyID,
                              time_block > window_start,
                              time_block <= window_end),
                       .(StudyID,
                         ema_time        = i.time_block,
                         functional_date = i.functional_date,
                         hr_missing,
                         sleep_missing),
                       nomatch = 0]
  
  window_matched <- matched[,
                            .(hr_missing_win   = mean(hr_missing,   na.rm = TRUE),
                              sleep_missing_win = as.integer(mean(sleep_missing, na.rm = TRUE) > 0)),
                            by = .(StudyID, ema_time, functional_date)]
  
  all_windows <- windows_dt[, .(StudyID, ema_time = time_block, functional_date)]
  window_missingness <- merge(all_windows, window_matched,
                              by = c("StudyID","ema_time","functional_date"),
                              all.x = TRUE)
  window_missingness[is.na(hr_missing_win),    hr_missing_win    := 1]
  window_missingness[is.na(sleep_missing_win), sleep_missing_win := 1L]
  window_missingness <- as.data.frame(window_missingness)
  
  window_eligibility <- window_missingness %>%
    group_by(StudyID) %>%
    summarise(
      n_ema_total        = n(),
      n_hr_valid_windows = sum(hr_missing_win <= 0.5, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      window_missingness %>%
        group_by(StudyID, functional_date) %>%
        summarise(sleep_missing_day = as.integer(any(sleep_missing_win > 0)),
                  .groups = "drop") %>%
        group_by(StudyID) %>%
        summarise(n_sleep_valid_days = sum(!sleep_missing_day), .groups = "drop"),
      by = "StudyID"
    )
  
  write_csv(window_eligibility, window_missingness_csv)
  cat("Saved: window_missingness_summary.csv\n")
}

cat("Window-based missingness summary:\n")
print(summary(window_eligibility[, c("n_ema_total","n_hr_valid_windows","n_sleep_valid_days")]))

eligible <- base_data %>%
  group_by(StudyID) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(window_eligibility, by = "StudyID") %>%
  filter(
    n_hr_valid_windows  >= hr_window_min,
    n_sleep_valid_days  >= sleep_day_min
  ) %>%
  arrange(desc(count)) %>%
  pull(StudyID)

cat("Eligible participants:", length(eligible),
    sprintf("(>= %d HR valid windows, >= %g sleep valid days)\n",
            hr_window_min, sleep_day_min))

add_opioid_lag1 <- function(ema_df, opioid_daily, prev_opioid = NULL, fill_na = NULL) {
  daily_lagged <- opioid_daily %>%
    arrange(functional_date) %>%
    mutate(opioid_num_lag1 = lag(opioid_num, 1))
  result <- ema_df %>%
    left_join(daily_lagged %>% select(functional_date, opioid_num_lag1),
              by = "functional_date")
  if (!is.null(prev_opioid)) {
    first_date <- min(ema_df$functional_date, na.rm = TRUE)
    result$opioid_num_lag1[result$functional_date == first_date &
                             is.na(result$opioid_num_lag1)] <- prev_opioid
  }
  if (!is.null(fill_na))
    result$opioid_num_lag1[is.na(result$opioid_num_lag1)] <- fill_na
  result
}

impute_fold <- function(train_5min, newdata_5min, m = m_imp, seed = 126) {
  mean_fb <- function(train_vec, new_vec) {
    mu <- mean(train_vec, na.rm = TRUE)
    if (is.nan(mu) || is.na(mu)) mu <- 0
    ifelse(is.na(new_vec), mu, new_vec)
  }
  vars_5min <- c("hr","steps")
  imp_cols  <- intersect(c("day","hour","minute","is_weekend","hr","steps"),
                         names(train_5min))
  d_train <- train_5min %>%
    select(all_of(imp_cols)) %>%
    mutate(across(where(is.logical), as.numeric))
  d_new <- newdata_5min %>%
    select(all_of(imp_cols)) %>%
    mutate(across(where(is.logical), as.numeric))
  
  imp_5min_tr <- NULL
  has_obs <- sapply(vars_5min[vars_5min %in% names(d_train)],
                    function(v) sum(!is.na(d_train[[v]])) > 1)
  if (length(has_obs) > 0 && all(has_obs)) {
    tryCatch({
      ini   <- mice(d_train, maxit = 0, printFlag = FALSE)
      pmeth <- ini$method
      pmat  <- ini$predictorMatrix
      n_days <- length(unique(d_train$day))
      for (v in vars_5min)
        if (v %in% names(pmeth) && pmeth[v] != "") {
          if (n_days >= 5 && "day" %in% names(d_train)) {
            pmeth[v]       <- "2l.pmm"
            pmat[v, "day"] <- -2
          } else {
            pmeth[v] <- "pmm"
          }
        }
      imp_5min_tr <- withCallingHandlers(
        mice(d_train, method=pmeth, predictorMatrix=pmat,
             m=m, maxit=5, ridge=0.1, seed=seed, printFlag=FALSE),
        warning = function(w) {
          if (grepl("longer object length", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    }, error=function(e) {
      message("  [2l.pmm failed, falling back to pmm] ", e$message)
      if (any(pmeth == "2l.pmm")) {
        tryCatch({
          pmeth_pmm <- pmeth
          pmeth_pmm[pmeth_pmm == "2l.pmm"] <- "pmm"
          imp_5min_tr <<- mice(d_train, method=pmeth_pmm,
                               predictorMatrix=ini$predictorMatrix,
                               m=m, maxit=5, ridge=0.1, seed=seed, printFlag=FALSE)
          message("  [pmm fallback succeeded for hr/steps]")
        }, error=function(e2) {
          message("  [pmm fallback also failed, will use mean] ", e2$message)
        })
      }
    })
  }
  if (is.null(imp_5min_tr) && !(length(has_obs) > 0 && all(has_obs)))
    message("  [hr/steps MICE skipped — insufficient observed values, using mean]")
  if (is.null(imp_5min_tr) && (length(has_obs) > 0 && all(has_obs)))
    message("  [hr/steps MICE failed entirely — using mean fill]")
  
  daily_tr <- train_5min %>%
    group_by(day) %>%
    summarise(
      totalMinutesAsleep = first_or_na(totalMinutesAsleep),
      opioid_num = if ("opioid_num" %in% names(train_5min)) first_or_na(opioid_num) else NA_real_,
      is_weekend = first(is_weekend),
      hr_mean    = mean(hr,    na.rm=TRUE),
      steps_sum  = sum(steps,  na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.logical), as.numeric),
           hr_mean   = ifelse(is.nan(hr_mean),   NA, hr_mean),
           steps_sum = ifelse(is.nan(steps_sum), NA, steps_sum))
  
  run_daily_mice <- function(target_col, label) {
    if (sum(!is.na(daily_tr[[target_col]])) <= 1) return(NULL)
    tryCatch(
      mice(daily_tr %>% select(day, is_weekend, hr_mean, steps_sum,
                               all_of(target_col)),
           method="pmm", m=m, maxit=5, seed=seed, printFlag=FALSE),
      error = function(e) { message("  mice ", label, " fold failed: ", e$message); NULL })
  }
  imp_sl_tr <- run_daily_mice("totalMinutesAsleep", "sleep")
  if (is.null(imp_sl_tr)) message("  [sleep MICE skipped or failed — using mean fill]")
  imp_op_tr <- run_daily_mice("opioid_num",         "opioid")
  if (is.null(imp_op_tr)) message("  [opioid MICE skipped or failed — using mean fill]")
  
  daily_nd <- newdata_5min %>%
    group_by(day) %>%
    summarise(
      totalMinutesAsleep = first_or_na(totalMinutesAsleep),
      opioid_num = if ("opioid_num" %in% names(newdata_5min)) first_or_na(opioid_num) else NA_real_,
      .groups = "drop"
    )
  daily_nd_full <- newdata_5min %>%
    select(day, functional_date) %>% distinct() %>%
    right_join(daily_nd, by = "day") %>%
    left_join(daily_tr %>% select(day, is_weekend, hr_mean, steps_sum), by = "day") %>%
    mutate(
      is_weekend = ifelse(is.na(is_weekend),
                          as.numeric(weekdays(functional_date) %in% c("Saturday","Sunday")),
                          is_weekend),
      hr_mean   = ifelse(is.na(hr_mean),   mean(daily_tr$hr_mean,   na.rm=TRUE), hr_mean),
      steps_sum = ifelse(is.na(steps_sum), mean(daily_tr$steps_sum, na.rm=TRUE), steps_sum)
    )
  
  mu_hr    <- mean(train_5min$hr,    na.rm=TRUE); if (is.nan(mu_hr)    || is.na(mu_hr))    mu_hr    <- 0
  mu_steps <- mean(train_5min$steps, na.rm=TRUE); if (is.nan(mu_steps) || is.na(mu_steps)) mu_steps <- 0
  mu_sl    <- mean(daily_tr$totalMinutesAsleep, na.rm=TRUE); if (is.nan(mu_sl) || is.na(mu_sl)) mu_sl <- 0
  mu_op    <- mean(daily_tr$opioid_num,         na.rm=TRUE); if (is.nan(mu_op) || is.na(mu_op)) mu_op <- 0
  
  result <- vector("list", m)
  for (i in seq_len(m)) {
    tr_i <- train_5min
    if (!is.null(imp_5min_tr)) {
      comp5 <- as.data.frame(complete(imp_5min_tr, i))
      for (v in vars_5min) if (v %in% names(comp5)) tr_i[[v]] <- comp5[[v]]
      n_na_hr    <- sum(is.na(tr_i$hr))
      n_na_steps <- sum(is.na(tr_i$steps))
      tr_i$hr    <- ifelse(is.na(tr_i$hr),    mu_hr,    tr_i$hr)
      tr_i$steps <- ifelse(is.na(tr_i$steps), mu_steps, tr_i$steps)
      if (n_na_hr    > 0) message(sprintf("  [post-MICE mean fill: %d hr NAs in imp %d]",    n_na_hr,    i))
      if (n_na_steps > 0) message(sprintf("  [post-MICE mean fill: %d steps NAs in imp %d]", n_na_steps, i))
    } else { for (v in vars_5min) tr_i[[v]] <- mean_fb(train_5min[[v]], tr_i[[v]]) }
    
    if (!is.null(imp_sl_tr)) {
      sl <- as.data.frame(complete(imp_sl_tr,i)) %>%
        select(day,totalMinutesAsleep) %>% rename(sl_imp=totalMinutesAsleep)
      tr_i <- tr_i %>% left_join(sl,by="day") %>%
        mutate(totalMinutesAsleep=sl_imp) %>% select(-sl_imp)
      tr_i$totalMinutesAsleep <- ifelse(is.na(tr_i$totalMinutesAsleep), mu_sl, tr_i$totalMinutesAsleep)
    } else {
      tr_i$totalMinutesAsleep <- ifelse(is.na(tr_i$totalMinutesAsleep), mu_sl, tr_i$totalMinutesAsleep)
    }
    if (!is.null(imp_op_tr)) {
      op <- as.data.frame(complete(imp_op_tr,i)) %>%
        select(day,opioid_num) %>% rename(op_imp=opioid_num)
      tr_i <- tr_i %>% left_join(op,by="day") %>%
        mutate(opioid_num=op_imp) %>% select(-op_imp)
      if ("opioid_num" %in% names(tr_i))
        tr_i$opioid_num <- ifelse(is.na(tr_i$opioid_num), mu_op, tr_i$opioid_num)
    } else if ("opioid_num" %in% names(tr_i)) {
      tr_i$opioid_num <- ifelse(is.na(tr_i$opioid_num), mu_op, tr_i$opioid_num)
    }
    
    nd_i <- newdata_5min
    if (!is.null(imp_5min_tr)) {
      tryCatch({
        reused <- mice.reuse(imp_5min_tr, d_new, maxit=3, printFlag=FALSE)
        comp_nd <- as.data.frame(reused[[i]])
        for (v in vars_5min) if (v %in% names(comp_nd)) nd_i[[v]] <- comp_nd[[v]]
        nd_i$hr    <- ifelse(is.na(nd_i$hr),    mu_hr,    nd_i$hr)
        nd_i$steps <- ifelse(is.na(nd_i$steps), mu_steps, nd_i$steps)
      }, error=function(e) {
        message("  [mice.reuse hr/steps failed, using mean] ", e$message)
        for (v in vars_5min) nd_i[[v]] <<- mean_fb(train_5min[[v]], nd_i[[v]]) })
    } else { for (v in vars_5min) nd_i[[v]] <- mean_fb(train_5min[[v]], nd_i[[v]]) }
    
    if (!is.null(imp_sl_tr)) {
      tryCatch({
        reused_sl <- mice.reuse(imp_sl_tr,
                                daily_nd_full %>% select(day,is_weekend,hr_mean,steps_sum,totalMinutesAsleep),
                                maxit=3, printFlag=FALSE)
        sl_nd <- as.data.frame(reused_sl[[i]]) %>%
          select(day,totalMinutesAsleep) %>% rename(sl_imp=totalMinutesAsleep)
        nd_i <- nd_i %>% left_join(sl_nd,by="day") %>%
          mutate(totalMinutesAsleep=ifelse(is.na(sl_imp),totalMinutesAsleep,sl_imp)) %>%
          select(-sl_imp)
        nd_i$totalMinutesAsleep <- ifelse(is.na(nd_i$totalMinutesAsleep), mu_sl, nd_i$totalMinutesAsleep)
      }, error=function(e) {
        message("  [mice.reuse sleep failed, using mean] ", e$message)
        nd_i$totalMinutesAsleep <<- ifelse(is.na(nd_i$totalMinutesAsleep), mu_sl, nd_i$totalMinutesAsleep) })
    }
    if (!is.null(imp_op_tr)) {
      tryCatch({
        reused_op <- mice.reuse(imp_op_tr,
                                daily_nd_full %>% select(day,is_weekend,hr_mean,steps_sum,opioid_num),
                                maxit=3, printFlag=FALSE)
        op_nd <- as.data.frame(reused_op[[i]]) %>%
          select(day,opioid_num) %>% rename(op_imp=opioid_num)
        nd_i <- nd_i %>% left_join(op_nd,by="day") %>%
          mutate(opioid_num=ifelse(is.na(op_imp),opioid_num,op_imp)) %>%
          select(-op_imp)
        if ("opioid_num" %in% names(nd_i))
          nd_i$opioid_num <- ifelse(is.na(nd_i$opioid_num), mu_op, nd_i$opioid_num)
      }, error=function(e) {
        message("  [mice.reuse opioid failed, using mean] ", e$message)
        if ("opioid_num" %in% names(nd_i))
          nd_i$opioid_num <<- ifelse(is.na(nd_i$opioid_num), mu_op, nd_i$opioid_num)
      })
    }
    result[[i]] <- list(train=tr_i, newdata=nd_i, mu_op=mu_op)
  }
  # return mice objects alongside completed data so impute_holdout can reuse them
  attr(result, "imp_5min_tr") <- imp_5min_tr
  attr(result, "imp_sl_tr")   <- imp_sl_tr
  attr(result, "imp_op_tr")   <- imp_op_tr
  attr(result, "mu_hr")       <- mu_hr
  attr(result, "mu_steps")    <- mu_steps
  attr(result, "mu_sl")       <- mu_sl
  attr(result, "mu_op")       <- mu_op
  attr(result, "d_train")     <- d_train
  result
}

impute_holdout <- function(last_fold_imp, nd5_holdout, m = m_imp) {
  # reuse trained mice objects from the last fold instead of re-running MICE
  imp_5min_tr <- attr(last_fold_imp, "imp_5min_tr")
  imp_sl_tr   <- attr(last_fold_imp, "imp_sl_tr")
  imp_op_tr   <- attr(last_fold_imp, "imp_op_tr")
  mu_hr       <- attr(last_fold_imp, "mu_hr")
  mu_steps    <- attr(last_fold_imp, "mu_steps")
  mu_sl       <- attr(last_fold_imp, "mu_sl")
  mu_op       <- attr(last_fold_imp, "mu_op")
  d_train     <- attr(last_fold_imp, "d_train")
  
  vars_5min <- c("hr","steps")
  mean_fb <- function(train_vec, new_vec) {
    mu <- mean(train_vec, na.rm = TRUE)
    if (is.nan(mu) || is.na(mu)) mu <- 0
    ifelse(is.na(new_vec), mu, new_vec)
  }
  
  d_new <- nd5_holdout %>%
    select(any_of(names(d_train))) %>%
    mutate(across(where(is.logical), as.numeric))
  
  daily_nd <- nd5_holdout %>%
    group_by(day) %>%
    summarise(
      totalMinutesAsleep = first_or_na(totalMinutesAsleep),
      opioid_num = if ("opioid_num" %in% names(nd5_holdout)) first_or_na(opioid_num) else NA_real_,
      .groups = "drop"
    )
  # use training means for hr_mean/steps_sum fill in daily_nd
  daily_tr_summary <- d_train %>%
    group_by(day) %>%
    summarise(hr_mean = mean(hr, na.rm=TRUE), steps_sum = sum(steps, na.rm=TRUE),
              is_weekend = first(is_weekend), .groups="drop")
  daily_nd_full <- nd5_holdout %>%
    select(day, functional_date) %>% distinct() %>%
    right_join(daily_nd, by = "day") %>%
    left_join(daily_tr_summary %>% select(day, is_weekend, hr_mean, steps_sum), by = "day") %>%
    mutate(
      is_weekend = ifelse(is.na(is_weekend),
                          as.numeric(weekdays(functional_date) %in% c("Saturday","Sunday")),
                          is_weekend),
      hr_mean   = ifelse(is.na(hr_mean),   mean(daily_tr_summary$hr_mean,   na.rm=TRUE), hr_mean),
      steps_sum = ifelse(is.na(steps_sum), mean(daily_tr_summary$steps_sum, na.rm=TRUE), steps_sum)
    )
  
  result <- vector("list", m)
  for (i in seq_len(m)) {
    nd_i <- nd5_holdout
    
    if (!is.null(imp_5min_tr)) {
      tryCatch({
        reused <- mice.reuse(imp_5min_tr, d_new, maxit=3, printFlag=FALSE)
        comp_nd <- as.data.frame(reused[[i]])
        for (v in vars_5min) if (v %in% names(comp_nd)) nd_i[[v]] <- comp_nd[[v]]
        nd_i$hr    <- ifelse(is.na(nd_i$hr),    mu_hr,    nd_i$hr)
        nd_i$steps <- ifelse(is.na(nd_i$steps), mu_steps, nd_i$steps)
      }, error = function(e) {
        message("  [holdout mice.reuse hr/steps failed, using mean] ", e$message)
        nd_i$hr    <<- ifelse(is.na(nd_i$hr),    mu_hr,    nd_i$hr)
        nd_i$steps <<- ifelse(is.na(nd_i$steps), mu_steps, nd_i$steps)
      })
    } else {
      nd_i$hr    <- ifelse(is.na(nd_i$hr),    mu_hr,    nd_i$hr)
      nd_i$steps <- ifelse(is.na(nd_i$steps), mu_steps, nd_i$steps)
    }
    
    if (!is.null(imp_sl_tr)) {
      tryCatch({
        reused_sl <- mice.reuse(imp_sl_tr,
                                daily_nd_full %>% select(day,is_weekend,hr_mean,steps_sum,totalMinutesAsleep),
                                maxit=3, printFlag=FALSE)
        sl_nd <- as.data.frame(reused_sl[[i]]) %>%
          select(day,totalMinutesAsleep) %>% rename(sl_imp=totalMinutesAsleep)
        nd_i <- nd_i %>% left_join(sl_nd, by="day") %>%
          mutate(totalMinutesAsleep=ifelse(is.na(sl_imp),totalMinutesAsleep,sl_imp)) %>%
          select(-sl_imp)
        nd_i$totalMinutesAsleep <- ifelse(is.na(nd_i$totalMinutesAsleep), mu_sl, nd_i$totalMinutesAsleep)
      }, error = function(e) {
        message("  [holdout mice.reuse sleep failed, using mean] ", e$message)
        nd_i$totalMinutesAsleep <<- ifelse(is.na(nd_i$totalMinutesAsleep), mu_sl, nd_i$totalMinutesAsleep)
      })
    }
    
    if (!is.null(imp_op_tr)) {
      tryCatch({
        reused_op <- mice.reuse(imp_op_tr,
                                daily_nd_full %>% select(day,is_weekend,hr_mean,steps_sum,opioid_num),
                                maxit=3, printFlag=FALSE)
        op_nd <- as.data.frame(reused_op[[i]]) %>%
          select(day,opioid_num) %>% rename(op_imp=opioid_num)
        nd_i <- nd_i %>% left_join(op_nd, by="day") %>%
          mutate(opioid_num=ifelse(is.na(op_imp),opioid_num,op_imp)) %>%
          select(-op_imp)
        if ("opioid_num" %in% names(nd_i))
          nd_i$opioid_num <- ifelse(is.na(nd_i$opioid_num), mu_op, nd_i$opioid_num)
      }, error = function(e) {
        message("  [holdout mice.reuse opioid failed, using mean] ", e$message)
        if ("opioid_num" %in% names(nd_i))
          nd_i$opioid_num <<- ifelse(is.na(nd_i$opioid_num), mu_op, nd_i$opioid_num)
      })
    }
    
    result[[i]] <- list(
      train   = last_fold_imp[[i]]$train,   # pass through unchanged
      newdata = nd_i,
      mu_op   = mu_op
    )
  }
  result
}

build_fold_features <- function(pid, ema_train, ema_val, ema_test,
                                fold_imp, m = m_imp) {
  get_feats <- function(ema_ts, sensor_5min)
    bind_rows(lapply(ema_ts, function(t) compute_ema_features(pid, t, sensor_5min)))
  
  get_opioid_daily <- function(df5) {
    if ("opioid_num" %in% names(df5))
      df5 %>% group_by(functional_date) %>%
      summarise(opioid_num = first_or_na(opioid_num), .groups = "drop")
    else
      df5 %>% group_by(functional_date) %>%
      summarise(opioid_num = NA_real_, .groups = "drop")
  }
  
  val_ts  <- ema_val$time_block
  test_ts <- ema_test$time_block
  nd_ts   <- c(val_ts, test_ts)
  
  merge_split <- function(ema_split, feat_df, op_daily, prev_op, fill_op) {
    ema_split %>%
      left_join(feat_df %>% select(-StudyID), by="time_block") %>%
      add_opioid_lag1(op_daily, prev_op, fill_na = fill_op) %>%
      arrange(time_block)
  }
  
  lapply(seq_len(m), function(ii) {
    tr5   <- fold_imp[[ii]]$train
    nd5   <- fold_imp[[ii]]$newdata
    mu_op <- fold_imp[[ii]]$mu_op
    
    feat_tr <- get_feats(ema_train$time_block, tr5)
    
    feat_nd  <- get_feats(nd_ts, nd5)
    feat_val <- feat_nd %>% filter(time_block %in% val_ts)
    feat_tst <- feat_nd %>% filter(time_block %in% test_ts)
    
    op_tr <- get_opioid_daily(tr5)
    op_nd <- get_opioid_daily(nd5)
    
    first_val_date <- min(ema_val$functional_date,  na.rm = TRUE)
    last_tr_op     <- op_tr %>%
      filter(functional_date == first_val_date - 1) %>%
      pull(opioid_num) %>% first_or_na()
    
    first_test_date <- min(ema_test$functional_date, na.rm = TRUE)
    last_val_op     <- op_nd %>%
      filter(functional_date == first_test_date - 1) %>%
      pull(opioid_num) %>% first_or_na()
    
    list(
      train = merge_split(ema_train, feat_tr,  op_tr, NULL,        mu_op),
      val   = merge_split(ema_val,   feat_val, op_nd, last_tr_op,  mu_op),
      test  = merge_split(ema_test,  feat_tst, op_nd, last_val_op, mu_op)
    )
  })
}

for (i in seq_along(eligible)) {
  pid          <- eligible[i]
  pid_fold_dir <- file.path(fold_cache_dir, pid)
  pid_base     <- base_data %>% filter(StudyID == pid) %>% arrange(time_block)
  
  ema_complete_all <- pid_base %>%
    filter(!is.na(pain_increasing),
           !is.na(overall_pain_lag2),
           !is.na(catastrophize_lag1), !is.na(catastrophize_lag2),
           !is.na(depress_lag1),       !is.na(depress_lag2),
           !is.na(interference_lag1),  !is.na(interference_lag2),
           !is.na(ema_missing_lag1))
  n_complete       <- nrow(ema_complete_all)
  ema_obs_complete <- ema_complete_all %>% slice(1:(n_complete - 5))
  ema_holdout      <- ema_complete_all %>% slice((n_complete - 4):n_complete)
  
  tmp_folds <- create_rolling_folds(ema_obs_complete, min_train=20, val_size=10,
                                    test_size=1, expanding=TRUE)
  if (length(tmp_folds) == 0) {
    cat(sprintf("[%d/%d] %s — no valid folds, skipping\n",i,length(eligible),pid))
    next
  }
  
  n_folds           <- length(tmp_folds)
  fold_rds_files_5b <- file.path(pid_fold_dir, sprintf("fold_%04d.rds", seq_len(n_folds)))
  holdout_rds_file  <- file.path(pid_fold_dir, "holdout.rds")
  
  if (all(file.exists(fold_rds_files_5b)) && file.exists(holdout_rds_file)) {
    cat(sprintf("[%d/%d] %s — fold RDS cached (%d folds)\n",
                i,length(eligible),pid,n_folds))
    next
  }
  
  # cache check failed - delete existing folder to remove orphaned files
  if (dir.exists(pid_fold_dir)) {
    unlink(pid_fold_dir, recursive = TRUE)
    cat(sprintf("  Deleted stale fold cache for %s\n", pid))
  }
  
  cat(sprintf("[%d/%d] %s — computing fold imputation + features (%d folds)\n",
              i,length(eligible),pid,n_folds))
  
  current_pid <- pid   # capture for use inside error handler
  
  tryCatch({
    if (!dir.exists(pid_fold_dir)) dir.create(pid_fold_dir, recursive=TRUE)
    
    pid_5min <- df_expanded %>%
      filter(StudyID == pid) %>%
      mutate(across(where(is.logical), as.numeric)) %>%
      arrange(time_block)
    
    for (fi in seq_along(tmp_folds)) {
      if (file.exists(fold_rds_files_5b[[fi]])) next
      fold <- tmp_folds[[fi]]
      
      ema_tr <- ema_obs_complete[fold$train,      ]
      ema_vl <- ema_obs_complete[fold$validation, ]
      ema_ts <- ema_obs_complete[fold$test,       ]
      
      train_end <- max(ema_tr$time_block)
      test_end  <- max(ema_ts$time_block)
      tr5  <- pid_5min %>% filter(time_block <= train_end)
      nd5  <- pid_5min %>% filter(time_block >  train_end & time_block <= test_end)
      
      if (nrow(tr5) < 2 || nrow(nd5) < 1) {
        cat(sprintf("  Fold %02d: insufficient 5-min rows, skipping\n", fi)); next }
      
      ctx_imp <- sprintf("[pid=%s] [fold=%04d] [impute_fold]", pid, fi)
      fold_imp <- withCallingHandlers(
        impute_fold(tr5, nd5, m=m_imp, seed=126+fi),
        warning = function(w) {
          if (!grepl("longer object length", conditionMessage(w)))
            write_log(log_imputation, "WARN", ctx_imp, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          write_log(log_imputation, "MSG", ctx_imp, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      ctx_feat <- sprintf("[pid=%s] [fold=%04d] [build_fold_features]", pid, fi)
      fold_features <- withCallingHandlers(
        build_fold_features(pid, ema_tr, ema_vl, ema_ts, fold_imp, m=m_imp),
        warning = function(w) {
          if (!grepl("longer object length", conditionMessage(w)))
            write_log(log_imputation, "WARN", ctx_feat, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          write_log(log_imputation, "MSG", ctx_feat, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      saveRDS(fold_features, fold_rds_files_5b[[fi]])
      cat(sprintf("  Saved fold %02d/%02d\n", fi, n_folds))
      
      # keep reference to last fold's imputation and training cutoff
      last_fold_imp       <- fold_imp
      last_fold_train_end <- train_end
    }
    
    holdout_end <- max(ema_holdout$time_block)
    nd5_holdout <- pid_5min %>%
      filter(time_block > last_fold_train_end & time_block <= holdout_end)
    
    if (nrow(nd5_holdout) >= 1) {
      ctx_ho <- sprintf("[pid=%s] [holdout] [impute_holdout]", pid)
      
      holdout_imp <- withCallingHandlers(
        impute_holdout(last_fold_imp, nd5_holdout, m=m_imp),
        warning = function(w) {
          if (!grepl("longer object length", conditionMessage(w)))
            write_log(log_imputation, "WARN", ctx_ho, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          write_log(log_imputation, "MSG", ctx_ho, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      ema_dummy_tr <- ema_obs_complete[nrow(ema_obs_complete), ]
      ema_dummy_vl <- ema_obs_complete[nrow(ema_obs_complete), ]
      holdout_features <- withCallingHandlers(
        build_fold_features(pid, ema_dummy_tr, ema_dummy_vl,
                            ema_holdout, holdout_imp, m=m_imp),
        warning = function(w) {
          if (!grepl("longer object length", conditionMessage(w)))
            write_log(log_imputation, "WARN", ctx_ho, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          write_log(log_imputation, "MSG", ctx_ho, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      
      holdout_rds <- lapply(holdout_features, function(x) x$test)
      saveRDS(holdout_rds, holdout_rds_file)
      cat("  Saved holdout with Fitbit + EMA features\n")
    } else {
      saveRDS(ema_holdout, holdout_rds_file)
      cat("  Saved holdout (raw EMA only — no sensor data for holdout window)\n")
    }
    
    gc()
    cat(sprintf("  Done: %s\n", pid))
  }, error = function(e) {
    msg <- sprintf("  Error processing %s in Section 5b: %s", current_pid, e$message)
    message(msg)
    write_log(log_imputation, "ERROR",
              sprintf("[pid=%s]", current_pid), e$message)
  })
}

cat("\n=== Fold-level feature pre-computation complete ===\n")


###############################################################################
# 6 - MAIN LOOP: ONE ANALYSIS AT A TIME, ONE PARTICIPANT AT A TIME
###############################################################################

write_log(log_modeling, "INFO", "[startup]", "Log for errors, warnings, and messages from models")

for (analysis in analyses) {
  an          <- analysis$name
  pvars       <- analysis$vars
  patient_dir <- file.path(results_dir, an, "PatientResults")
  summary_dir <- file.path(analysis_output_dir, an, "Summary")
  
  cat("\n\n========================================\n")
  cat("Analysis:", an, "\n")
  cat("========================================\n")
  
  summary_csv <- file.path(summary_dir, paste0(an, "_participant_summary.csv"))
  
  # check results_dir for existing participant results
  has_results <- sapply(eligible, function(pid) {
    file.exists(file.path(patient_dir, pid,
                          paste0(pid, "_participant_summary_row.rds")))
  })
  done_pids <- eligible[has_results]
  remaining <- eligible[!has_results]
  
  cat("Participants with existing results:", length(done_pids), "\n")
  cat("Remaining to process:", length(remaining), "\n")
  
  existing_summary <- bind_rows(lapply(done_pids, function(pid) {
    readRDS(file.path(patient_dir, pid,
                      paste0(pid, "_participant_summary_row.rds")))
  }))
  
  for (i in seq_along(remaining)) {
    pid <- remaining[i]
    cat("\n[", an, "] Processing", i, "/", length(remaining), ":", pid, "\n")
    set.seed(126 + which(eligible == pid))
    
    pid_fold_dir     <- file.path(fold_cache_dir, pid)
    fold_rds_files   <- sort(list.files(pid_fold_dir, pattern="^fold_.*\\.rds$",
                                        full.names = TRUE))
    holdout_rds_file <- file.path(pid_fold_dir, "holdout.rds")
    if (length(fold_rds_files) == 0) {
      cat("  No fold RDS files for", pid, "(skipped)\n")
      next
    }
    cat(sprintf("  Loading %d fold RDS files\n", length(fold_rds_files)))
    
    ctx_pid <- sprintf("[%s] [pid=%s]", an, pid)
    res <- withCallingHandlers(
      tryCatch(
        build_participant_model(pid, pvars, patient_dir,
                                analysis_name    = an,
                                fold_rds_files   = fold_rds_files,
                                holdout_rds_file = holdout_rds_file),
        error = function(e) {
          message("  Unhandled error for ", pid, ": ", e$message)
          write_log(log_modeling, "ERROR", ctx_pid, e$message)
          NULL
        }
      ),
      warning = function(w) {
        write_log(log_modeling, "WARN", ctx_pid, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        write_log(log_modeling, "MSG", ctx_pid, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    
    if (!is.null(res)) {
      existing_summary <- bind_rows(existing_summary, res)
      write_csv(existing_summary, summary_csv)
      cat("  Saved summary row for", pid,
          "| Ensemble AUC:", round(res$Ensemble_AUC, 3), "\n")
    } else {
      cat("  No result for", pid, "(skipped)\n")
    }
  }
  
  # write summary CSVs
  if (nrow(existing_summary) > 0)
    write_csv(existing_summary, summary_csv)
  
  cat("\n[", an, "] Complete.",
      nrow(existing_summary), "participants in summary.\n")
}

cat("\n=== Stage 2 complete. Run step3_analysis.R ===\n")
cat(sprintf("Imputation log: %s\n", log_imputation))
cat(sprintf("Modeling log:   %s\n", log_modeling))
