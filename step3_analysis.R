################################################################################
#
# Personalized Pain Prediction (PPP) - Step 3: Aggregation + Visualization
#
#   Reads the per-participant CSVs written by step2_imputation_and_modeling.R, assembles
#   summary tables, runs hypothesis tests, and produces all publication
#   plots. 
#
# Run after step2_imputation_and_modeling.R.
#
################################################################################


###############################################################################
# 0 - LIBRARIES
###############################################################################

pkgs <- c("readr","dplyr","tidyr","ggplot2","pROC","grid",
          "ggridges","patchwork","forcats")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


###############################################################################
# 1 - PATHS
###############################################################################

source("paths_local.R")
hr_window_min  <- 26   # must match the value used in step2_imputation_and_modeling.R
analysis_label <- paste0(hr_window_min, "ValidWindows")

analysis_output_dir <- file.path(output_dir, analysis_label)
results_dir         <- file.path(output_dir, "results")

analyses <- list(
  list(name = "A1_FitbitOnly",   label = "Fitbit + Context"),
  list(name = "A2_LagEMAOnly",   label = "Lag EMA Only"),
  list(name = "A3_Combined",     label = "Fitbit + Lag EMA + Context"),
  list(name = "A4_Lag1PainOnly", label = "Lag-1 Pain Only") 
)


###############################################################################
# 2 - UTILITY FUNCTIONS
###############################################################################

# Saves a plot as both PNG (for quick viewing) and TIFF (for publication
# submission), at the same publication-scaled dimensions and dpi = 300.
save_fig <- function(path_no_ext, plot, width, height, dpi = 300) {
  ggsave(paste0(path_no_ext, ".png"),  plot = plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(path_no_ext, ".tiff"), plot = plot, width = width, height = height, dpi = dpi,
         compression = "lzw")
}

format_decimal <- function(x, digits = 3) {
  xr <- round(x, digits)
  xc <- as.character(xr)
  xc <- ifelse(abs(xr - 0) < 1e-10, "0", xc)
  xc <- ifelse(abs(xr - 1) < 1e-10, "1", xc)
  xc
}

clean_variable_names <- function(var_names) {
  mapping <- c(
    "overall_pain_lag1"   = "Overall Pain Lag1",
    "overall_pain_lag2"   = "Overall Pain Lag2",
    "catastrophize_lag1"  = "Catastrophize Lag1",
    "catastrophize_lag2"  = "Catastrophize Lag2",
    "depress_lag1"        = "Depression Lag1",
    "depress_lag2"        = "Depression Lag2",
    "interference_lag1"   = "Interference Lag1",
    "interference_lag2"   = "Interference Lag2",
    "opioid_num_lag1"     = "Opioid Number Lag1",
    "opioid_num_lag2"     = "Opioid Number Lag2",
    "ema_missing"         = "EMA Missed",
    "averageHR"           = "Average HR",
    "varianceHR"          = "Variance HR",
    "wearTime"            = "Wear Time (mins)",
    "totalSteps"          = "Total Steps",
    "averageSteps"        = "Average Steps",
    "varianceSteps"       = "Variance Steps",
    "sedentaryTime"       = "Sedentary Time (mins)",
    "activeTime"          = "Active Time (mins)",
    "totalMinutesAsleep"  = "Total Minutes Asleep",
    "hr_missing"          = "Proportion Fitbit Missing",
    "steps_missing"       = "Proportion Steps Missing",
    "sleep_missing"       = "Proportion Sleep Missing",
    "weather_temp"        = "Ambient Temperature",
    "weather_humidity"    = "Humidity",
    "weather_pressure"    = "Barometric Pressure",
    "weather_precip"      = "Precipitation",
    "is_weekend"          = "Weekend",
    "hour"                = "Hour",
    "days_since_first_ema" = "Days Since First EMA"
  )
  ifelse(var_names %in% names(mapping), mapping[var_names], var_names)
}

calc_ci <- function(x) {
  n <- sum(!is.na(x)); if (n < 2) return(c(NA, NA))
  m  <- mean(x, na.rm = TRUE); se <- sd(x, na.rm = TRUE) / sqrt(n)
  tc <- qt(.975, n - 1); c(m - tc*se, m + tc*se)
}

fmt_est_sd_ci <- function(m, s, lo, hi)
  paste0(round(m, 3), " (", round(s, 3), ") [", round(lo, 3), ", ", round(hi, 3), "]")

get_en_col <- function(sm) {
  if ("LogisticRegression_AUC" %in% names(sm)) "LogisticRegression_AUC"
  else "ElasticNet_AUC"
}

get_auc_cols <- function(sm) {
  en_col <- get_en_col(sm)
  c("RandomForest_AUC", en_col, "GaussianProcess_AUC", "Ensemble_AUC")
}

get_en_label <- function(sm) {
  if ("LogisticRegression_AUC" %in% names(sm)) "Logistic Regression"
  else "Elastic Net"
}


###############################################################################
# 3 - LOAD ALL PER-PARTICIPANT RESULTS FROM DISK
###############################################################################

cat("=== Stage 3: Loading results from disk ===\n\n")

all_summaries  <- list()
all_sens_spec  <- list()
all_importance <- list()

for (a in analyses) {
  an          <- a$name
  summary_csv <- file.path(analysis_output_dir, an, "Summary",
                           paste0(an, "_participant_summary.csv"))
  patient_dir <- file.path(results_dir, an, "PatientResults")
  
  if (file.exists(summary_csv)) {
    sm <- read_csv(summary_csv, show_col_types = FALSE)
    all_summaries[[an]] <- sm
    cat(an, "— summary rows:", nrow(sm), "\n")
  } else {
    cat(an, "— no summary CSV found (stage 2 not run yet?)\n")
    next
  }
  
  ss_list  <- list()
  imp_list <- list()
  
  for (pid in sm$StudyID) {
    pdir <- file.path(patient_dir, pid)
    ss_file  <- file.path(pdir, paste0(pid, "_pooled_sens_spec.csv"))
    imp_file <- file.path(pdir, paste0(pid, "_feature_importance.csv"))
    if (file.exists(ss_file))
      ss_list[[pid]] <- read_csv(ss_file, show_col_types = FALSE)
    if (file.exists(imp_file))
      imp_list[[pid]] <- read_csv(imp_file, show_col_types = FALSE) %>%
      mutate(ParticipantID = pid)
  }
  
  if (length(ss_list)  > 0) all_sens_spec[[an]]  <- bind_rows(ss_list)
  if (length(imp_list) > 0) all_importance[[an]] <- bind_rows(imp_list)
}


###############################################################################
# 3b - DETECT PARTICIPANTS WITH DEGENERATE MODEL CIs (GLOBAL, ACROSS ANALYSES)
###############################################################################

cat("\n--- Detecting participants with degenerate CIs ---\n")

# Parse "[lo, hi]" strings into numeric, flag [NA,NA] / [0,1] / [1,1] as degenerate
parse_ci <- function(x) suppressWarnings(as.numeric(strsplit(gsub("\\[|\\]", "", x), ",\\s*")[[1]]))
is_bad   <- function(ci) (is.na(ci[1]) && is.na(ci[2])) || isTRUE(all.equal(ci, c(0,1))) || isTRUE(all.equal(ci, c(1,1)))

# Step 1: null out any model (RF, EN/LR, GP, Ensemble) whose own CI is degenerate
all_summaries <- lapply(all_summaries, function(sm) {
  pairs <- list(c("RandomForest_AUC","RF_CI"), 
                c(get_en_col(sm), if("LR_CI" %in% names(sm)) "LR_CI" else "EN_CI"),
                c("GaussianProcess_AUC","GP_CI"), 
                c("Ensemble_AUC","Ens_CI"))
  for (p in pairs) {
    bad <- sapply(sm[[p[2]]], function(x) is_bad(parse_ci(x)))
    sm[[p[1]]][bad] <- NA
    sm[[p[2]]][bad] <- "[NA, NA]"
  }
  sm
})

# Step 2: find participants with zero usable models in ANY config (global)
auc_cols_by_sm <- lapply(all_summaries, function(sm) c("RandomForest_AUC", get_en_col(sm), "GaussianProcess_AUC", "Ensemble_AUC"))
drop_ids_degenerate <- unique(unlist(lapply(names(all_summaries), function(an) {
  sm <- all_summaries[[an]]
  sm$StudyID[apply(sm[, auc_cols_by_sm[[an]]], 1, function(r) all(is.na(r)))]
})))

# Step 2b: find participants with only a single held-out observation in ANY
# config (global)
drop_ids_single_obs <- unique(unlist(lapply(names(all_summaries), function(an) {
  sm <- all_summaries[[an]]
  # remove participants with only one fold
  sm$StudyID[sm$NumObservations == 1]
})))

drop_ids <- unique(c(drop_ids_degenerate, drop_ids_single_obs))

# Step 3: drop those participants from every config
all_summaries <- lapply(all_summaries, function(sm) sm %>% filter(!StudyID %in% drop_ids))


###############################################################################
# 4 - AGGREGATE FEATURE IMPORTANCE PER ANALYSIS
###############################################################################

cat("\n--- Aggregating feature importance ---\n")

for (a in analyses) {
  an  <- a$name
  imp <- all_importance[[an]]
  if (is.null(imp) || nrow(imp) == 0) next
  
  sm          <- all_summaries[[an]]
  summary_dir <- file.path(analysis_output_dir, an, "Summary")
  auc_cols    <- get_auc_cols(sm)
  en_label    <- get_en_label(sm)
  model_order <- c("Random Forest", en_label, "Gaussian Process", "Ensemble")
  
  top_list <- list()
  for (pid in sm$StudyID) {
    pr  <- sm[sm$StudyID == pid, ]
    bm  <- model_order[which.max(unlist(pr[, auc_cols]))]
    pid_imp <- imp %>% filter(ParticipantID == pid, Model == bm)
    if (nrow(pid_imp) == 0) next
    tf <- pid_imp %>% arrange(desc(Mean_Importance)) %>% head(3) %>%
      mutate(Variable_Clean = clean_variable_names(Variable)) %>%
      pull(Variable_Clean)
    if (length(tf) > 0)
      top_list[[pid]] <- data.frame(ParticipantID = pid, BestModel = bm,
                                    Feature = tf, Rank = seq_along(tf))
  }
  
  if (length(top_list) > 0) {
    all_tf <- bind_rows(top_list)
    n_participants <- length(unique(all_tf$ParticipantID))
    
    fc <- all_tf %>%
      count(Feature) %>%
      mutate(Pct = round(n / n_participants * 100, 1)) %>%
      arrange(desc(Pct))
    
    write_csv(all_tf, file.path(summary_dir, paste0(an,"_individual_top_features.csv")))
    write_csv(fc,     file.path(summary_dir, paste0(an,"_top_features_summary.csv")))
    all_importance[[an]] <- fc
    cat(an, "— top features saved\n")
  }
}


###############################################################################
# 5 - VISUALISATION HELPERS
###############################################################################

colors <- c("Random Forest"      = "#2A363B",
            "Elastic Net"         = "#83AF9B",
            "Logistic Regression" = "#83AF9B",
            "Gaussian Process"    = "#E6C76E",
            "Ensemble"            = "#FE4365")

.build_three_plots <- function(sm, label) {
  
  en_col   <- get_en_col(sm)
  en_label <- get_en_label(sm)
  
  # pivot the four AUC columns, handling either EN or LR column naming
  long <- sm %>%
    select(StudyID,
           RandomForest_AUC,
           !!sym(en_col),
           GaussianProcess_AUC,
           Ensemble_AUC) %>%
    rename(EN_AUC = !!sym(en_col)) %>%
    pivot_longer(-StudyID, names_to = "Model", values_to = "AUC") %>%
    mutate(Model = case_match(Model,
                              "RandomForest_AUC"    ~ "Random Forest",
                              "EN_AUC"              ~ en_label,
                              "GaussianProcess_AUC" ~ "Gaussian Process",
                              "Ensemble_AUC"        ~ "Ensemble",
                              .default = Model),
           Model = factor(Model, levels = names(colors)))
  
  above <- long %>% group_by(StudyID) %>%
    filter(any(AUC > 0.5)) %>% pull(StudyID) %>% unique()
  filt  <- long %>% filter(StudyID %in% above)
  
  auc_cols <- get_auc_cols(sm)
  best <- sm %>% rowwise() %>%
    mutate(
      Best_AUC   = max(!!!syms(auc_cols), na.rm = TRUE),
      Best_Model = {
        aucs <- c_across(all_of(auc_cols))
        if (all(is.na(aucs))) NA_character_
        else c("Random Forest", en_label, "Gaussian Process", "Ensemble")[which.max(aucs)]
      }) %>%
    ungroup() %>%
    filter(StudyID %in% above) %>%
    mutate(Best_Model = factor(Best_Model, levels = names(colors)))
  
  filt_ord <- filt %>% left_join(best %>% select(StudyID, Best_AUC), by = "StudyID")
  
  p_ridge <- ggplot(filt, aes(x = AUC, y = Model, fill = Model)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9, rel_min_height = 0.01,
                        aes(color = Model), quantile_lines = TRUE, jittered_points = TRUE) +
    scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
    scale_y_discrete(expand = expansion(mult = c(.1, .1))) +
    labs(x = paste0("Pooled ROC AUC (", label, ")"), y = "") +
    theme_minimal() +
    theme(axis.ticks = element_blank(), legend.position = "none",
          axis.text = element_text(size = 6), axis.title.x = element_text(size = 6)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25),
                       labels = format_decimal(seq(0,1,.25)))

  p_dot <- ggplot(filt_ord, aes(x = AUC, y = fct_reorder(StudyID, Best_AUC),
                                color = Model)) +
    geom_point(size = .2, alpha = 0.8) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = 0.5,  linetype = "dashed", alpha = 0.5, linewidth=0.25) +
    labs(x = paste0("Pooled ROC AUC (", label, ")"),
         y = "Participants (ranked by best performance)", color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_blank(), axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6), legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"), legend.key.spacing = unit(0.2, "cm"),
          legend.key.spacing.y = unit(0.01, "cm")) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25),
                       labels = format_decimal(seq(0,1,.25)))

  freq <- best %>% count(Best_Model) %>% mutate(Pct = round(n / sum(n) * 100, 1))
  p_freq <- ggplot(freq, aes(x = n, y = reorder(Best_Model, n), fill = Best_Model)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = colors) +
    labs(x = paste0("Number of Participants (", label, ")"), y = "") +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(),
          axis.text = element_text(size = 6), axis.title.x = element_text(size = 6),
          plot.margin = unit(c(.25,.8,.25,.5), "cm")) +
    geom_text(aes(label = n), hjust = 1.2, size = 5 / .pt, color = "white", fontface = "bold") +
    geom_text(aes(label = paste0("(", format_decimal(Pct), "%)")),
              hjust = -0.1, size = 5 / .pt, color = "black", fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, .2))) + coord_cartesian(clip = "off")
  
  list(p_dot = p_dot, p_freq = p_freq, p_ridge = p_ridge)
}


###############################################################################
# 6 - COMBINED MODEL PERFORMANCE PLOTS
###############################################################################

cat("\n--- Creating model performance plots ---\n")

af    <- Filter(function(a) a$name != "A4_Lag1PainOnly", analyses)
plots <- list()

for (i in seq_along(af)) {
  sm <- all_summaries[[af[[i]]$name]]
  if (!is.null(sm) && nrow(sm) > 0) {
    pp <- .build_three_plots(sm, af[[i]]$label)
    pp$p_dot <- pp$p_dot +
      annotation_custom(grid::textGrob(LETTERS[i], x = .02, y = .98, hjust = 0, vjust = 1,
                                       gp = grid::gpar(fontsize = 7, fontface = "bold")))
    plots[[af[[i]]$name]] <- pp
  }
}

if (length(plots) > 0) {
  p_dots <- wrap_plots(lapply(af, function(a) plots[[a$name]]$p_dot),
                       ncol = length(af)) +
    plot_layout(guides = "collect") &
    theme(text = element_text(family = "Arial"),
          legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          legend.key.spacing = unit(0.2, "cm"),
          legend.key.spacing.y = unit(0.01, "cm"))
  save_fig(file.path(analysis_output_dir, "combined_participant_plots"),
           plot = p_dots, width = 7, height = 3)

  p_freqs <- wrap_plots(lapply(af, function(a) plots[[a$name]]$p_freq),
                        ncol = length(af)) &
    theme(text = element_text(family = "Arial"),
          plot.margin = unit(c(.25,.9,.25,.1), "cm"))
  save_fig(file.path(analysis_output_dir, "combined_frequency_plots"),
           plot = p_freqs, width = 6, height = 2.5)
  cat("  Saved: combined_participant_plots.png/.tiff, combined_frequency_plots.png/.tiff\n")
}

a4    <- Filter(function(a) a$name == "A4_Lag1PainOnly", analyses)[[1]]
sm_a4 <- all_summaries[[a4$name]]

if (!is.null(sm_a4) && nrow(sm_a4) > 0) {
  pp  <- .build_three_plots(sm_a4, a4$label)
  # Drop the chance-line (xintercept = 0.5) baked into .build_three_plots()'s
  # p_dot
  pp$p_dot$layers <- Filter(function(l) !(inherits(l$geom, "GeomVline") &&
                                             isTRUE(all.equal(l$data$xintercept, 0.5))),
                             pp$p_dot$layers)
  # Reference line for the null-model comparison (mean AUC = 0.557)
  pp$p_dot <- pp$p_dot +
    geom_vline(xintercept = 0.557, linetype = "dashed", color = "grey40", alpha = 0.8, linewidth=0.25)
  p_a4 <- ((pp$p_dot + pp$p_freq + plot_layout(widths = c(2,1))) &
             theme(text = element_text(family = "Arial"))) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 7, face = "bold"),
          plot.tag.position = c(0, 1))
  save_fig(file.path(analysis_output_dir, "Lag1PainOnly_combined_grid"),
           plot = p_a4, width = 7, height = 3)
  cat("  Saved: Lag1PainOnly_combined_grid.png/.tiff\n")
}

###############################################################################
# 7 - COMBINED FEATURE IMPORTANCE PLOTS
###############################################################################

cat("\n--- Creating feature importance plots ---\n")

mc_colors  <- c("#2A363B","#83AF9B","#E6C76E","#FE4365")
plist      <- list()
af_imp     <- Filter(function(a) a$name != "A4_Lag1PainOnly", analyses)
letter_idx <- 0

for (i in seq_along(af_imp)) {
  an  <- af_imp[[i]]$name
  imp <- all_importance[[an]]
  if (!is.null(imp) && nrow(imp) > 0 && "Pct" %in% names(imp)) {
    letter_idx <- letter_idx + 1
    impd <- imp %>% arrange(desc(Pct)) %>%
      slice_head(n = 10) %>% mutate(Feature = factor(Feature, levels = rev(Feature)))
    write_csv(impd, file.path(analysis_output_dir, paste0(an, "_imp_data.csv")))
    plist[[an]] <- ggplot(impd, aes(x = Pct, y = Feature)) +
      geom_bar(stat = "identity", fill = mc_colors[i], alpha = 0.8) +
      labs(x = paste0("Percentage (", af_imp[[i]]$label, ")"), y = "") +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            plot.margin = unit(c(1,1,.5,.5), "cm")) +
      annotation_custom(grid::textGrob(LETTERS[letter_idx], x = -.08, y = 1.02,
                                       hjust = 0, vjust = 0,
                                       gp = grid::gpar(fontsize = 7, fontface = "bold"))) +
      coord_cartesian(clip = "off")
  }
}

if (length(plist) > 0) {
  p_imp_combined <- wrap_plots(plist, ncol = length(plist)) &
    theme(plot.margin = unit(c(.3,.25,.25,.25), "cm"))
  save_fig(file.path(analysis_output_dir, "combined_feature_importance"),
           plot = p_imp_combined, width = 8, height = 2)
  cat("  Saved: combined_feature_importance.png/.tiff\n")
}


###############################################################################
# 8 - HYPOTHESIS TESTING (paired t-tests on pooled AUCs)
###############################################################################

cat("\n--- Hypothesis testing ---\n")

comps <- list(
  list(n1="A1_FitbitOnly",   n2="A2_LagEMAOnly",    l1="A1", l2="A2"),
  list(n1="A1_FitbitOnly",   n2="A3_Combined",      l1="A1", l2="A3"),
  list(n1="A1_FitbitOnly",   n2="A4_Lag1PainOnly",  l1="A1", l2="A4"),
  list(n1="A2_LagEMAOnly",   n2="A3_Combined",      l1="A2", l2="A3"),
  list(n1="A2_LagEMAOnly",   n2="A4_Lag1PainOnly",  l1="A2", l2="A4"),
  list(n1="A3_Combined",     n2="A4_Lag1PainOnly",  l1="A3", l2="A4")
)

ht_results <- data.frame()
for (co in comps) {
  d1 <- all_summaries[[co$n1]]; d2 <- all_summaries[[co$n2]]
  if (is.null(d1) || is.null(d2)) next
  ids <- intersect(d1$StudyID, d2$StudyID)
  if (length(ids) <= 5) next
  d1m <- d1[d1$StudyID %in% ids, ] %>% arrange(StudyID)
  d2m <- d2[d2$StudyID %in% ids, ] %>% arrange(StudyID)
  # use get_auc_cols to handle EN vs LR column name per analysis
  b1  <- apply(d1m[, get_auc_cols(d1m)], 1, max, na.rm = TRUE)
  b2  <- apply(d2m[, get_auc_cols(d2m)], 1, max, na.rm = TRUE)
  tt  <- t.test(b1, b2, paired = TRUE)
  dif <- mean(b1 - b2); n <- length(b1); se <- sd(b1 - b2) / sqrt(n)
  tc  <- qt(.975, n - 1)
  ht_results <- bind_rows(ht_results, data.frame(
    Comparison   = paste(co$l1, "vs", co$l2),
    Mean_AUC_1   = round(mean(b1, na.rm = TRUE), 3),
    Mean_AUC_2   = round(mean(b2, na.rm = TRUE), 3),
    Difference   = round(dif, 3),
    CI_Lower     = round(dif - tc*se, 3),
    CI_Upper     = round(dif + tc*se, 3),
    t_statistic  = round(tt$statistic, 2),
    P_value_raw  = tt$p.value,
    P_value      = ifelse(tt$p.value < .001, "<0.001",
                          format_decimal(tt$p.value, 3))
  ))
}

if (nrow(ht_results) > 0) {
  ht_results$Diff_CI <- paste0(ht_results$Difference, " (",
                               ht_results$CI_Lower, ", ",
                               ht_results$CI_Upper, ")")
  ht_results$P_adj         <- p.adjust(ht_results$P_value_raw, method = "BH")
  ht_results$P_adj_display <- ifelse(ht_results$P_adj < .001, "<0.001",
                                     format_decimal(ht_results$P_adj, 3))
} else {
  ht_results$Diff_CI       <- character(0)
  ht_results$P_adj         <- numeric(0)
  ht_results$P_adj_display <- character(0)
}
write_csv(ht_results, file.path(analysis_output_dir, "paired_t_test_pooled_auc.csv"))
cat("  Saved: paired_t_test_pooled_auc.csv\n")


###############################################################################
# 9 - SENSITIVITY / SPECIFICITY SUMMARY TABLE
###############################################################################

cat("\n--- Sensitivity / specificity summary ---\n")

combined_ss <- data.frame()
for (a in analyses) {
  ss <- all_sens_spec[[a$name]]
  if (!is.null(ss) && nrow(ss) > 0) {
    ss$Analysis       <- a$name
    ss$Analysis_Label <- a$label
    combined_ss <- bind_rows(combined_ss, ss)
  }
}

if (nrow(combined_ss) > 0) {
  write_csv(combined_ss, file.path(analysis_output_dir, "sensitivity_specificity_all.csv"))
  
  best_ss <- data.frame()
  for (a in analyses) {
    ss <- combined_ss %>% filter(Analysis == a$name)
    sm <- all_summaries[[a$name]]
    if (is.null(sm) || nrow(ss) == 0) next
    auc_cols    <- get_auc_cols(sm)
    en_label    <- get_en_label(sm)
    model_order <- c("Random Forest", en_label, "Gaussian Process", "Ensemble")
    for (pid in unique(ss$StudyID)) {
      pr  <- sm[sm$StudyID == pid, ]
      if (nrow(pr) == 0) next
      bm  <- model_order[which.max(unlist(pr[, auc_cols]))]
      br  <- ss %>% filter(StudyID == pid, Model == bm)
      if (nrow(br) > 0) best_ss <- bind_rows(best_ss, br)
    }
  }
  
  tbl2 <- best_ss %>%
    group_by(Analysis_Label) %>%
    summarise(
      N           = n(),
      AUC         = fmt_est_sd_ci(mean(AUC,         na.rm=TRUE), sd(AUC,         na.rm=TRUE), calc_ci(AUC)[1],         calc_ci(AUC)[2]),
      Sensitivity = fmt_est_sd_ci(mean(Sensitivity, na.rm=TRUE), sd(Sensitivity, na.rm=TRUE), calc_ci(Sensitivity)[1], calc_ci(Sensitivity)[2]),
      Specificity = fmt_est_sd_ci(mean(Specificity, na.rm=TRUE), sd(Specificity, na.rm=TRUE), calc_ci(Specificity)[1], calc_ci(Specificity)[2]),
      PPV         = fmt_est_sd_ci(mean(PPV,         na.rm=TRUE), sd(PPV,         na.rm=TRUE), calc_ci(PPV)[1],         calc_ci(PPV)[2]),
      NPV         = fmt_est_sd_ci(mean(NPV,         na.rm=TRUE), sd(NPV,         na.rm=TRUE), calc_ci(NPV)[1],         calc_ci(NPV)[2]),
      F1          = fmt_est_sd_ci(mean(F1,          na.rm=TRUE), sd(F1,          na.rm=TRUE), calc_ci(F1)[1],          calc_ci(F1)[2]),
      Kappa       = fmt_est_sd_ci(mean(Kappa,       na.rm=TRUE), sd(Kappa,       na.rm=TRUE), calc_ci(Kappa)[1],       calc_ci(Kappa)[2]),
      .groups     = "drop"
    )
  colnames(tbl2)[1] <- "Feature Configuration"
  
  write_csv(tbl2,    file.path(analysis_output_dir, "Table2_performance_metrics_pooled.csv"))
  write_csv(best_ss, file.path(analysis_output_dir, "best_model_sens_spec_detail.csv"))
  cat("  Saved: Table2_performance_metrics_pooled.csv\n")
  print(tbl2 %>% select(`Feature Configuration`, N, AUC, Sensitivity, Specificity))
}


###############################################################################
# 10 - PROGRESS REPORT
###############################################################################

cat("\n=== Stage 3 complete ===\n\n")
cat("Participant counts per analysis:\n")
for (a in analyses) {
  sm <- all_summaries[[a$name]]
  cat(sprintf("  %-30s  %d participants\n", a$name,
              if (!is.null(sm)) nrow(sm) else 0))
}
cat("\nAll outputs written to:", analysis_output_dir, "\n")


###############################################################################
# 11 - A4 LAG-1 PAIN LR COEFFICIENT SIGN
###############################################################################
#
#   Reports the sign and magnitude of each participant's A4_Lag1PainOnly LR
#   coefficient on overall_pain_lag1, using the value already computed by
#   build_a4_result_row() in step2.
#
###############################################################################

cat("\n--- A4 lag-1 pain LR coefficient sign ---\n")

sm_a4_lr <- all_summaries[["A4_Lag1PainOnly"]]

if (is.null(sm_a4_lr) || !("Lag1Pain_LR_Estimate" %in% names(sm_a4_lr))) {
  cat("  A4_Lag1PainOnly summary missing or has no Lag1Pain_LR_Estimate column — skipping\n")
} else {
  lr_signs <- sm_a4_lr %>%
    filter(!is.na(Lag1Pain_LR_Estimate)) %>%
    transmute(
      StudyID  = StudyID,
      Estimate = Lag1Pain_LR_Estimate,
      SE       = Lag1Pain_LR_SE,
      CI_Lo    = Lag1Pain_LR_CI_Lo,
      CI_Hi    = Lag1Pain_LR_CI_Hi,
      P        = Lag1Pain_LR_P,
      Sign     = ifelse(Estimate > 0, "Positive", "Negative")
    )

  sign_summary <- lr_signs %>%
    count(Sign) %>%
    mutate(Pct = round(n / sum(n) * 100, 1))

  write_csv(lr_signs,     file.path(analysis_output_dir, "A4_Lag1PainOnly_LR_coefficient_signs.csv"))
  write_csv(sign_summary, file.path(analysis_output_dir, "A4_Lag1PainOnly_LR_coefficient_sign_summary.csv"))

  cat("  Saved: A4_Lag1PainOnly_LR_coefficient_signs.csv, A4_Lag1PainOnly_LR_coefficient_sign_summary.csv\n")
  cat("  N with a usable coefficient:", nrow(lr_signs), "/", nrow(sm_a4_lr), "\n")
  print(sign_summary)
}


###############################################################################
# 12 - A4 NULL MODEL COMPARISON (permutation test)
###############################################################################
#
#   Reads the per-participant null distributions written by the A4_NullModel
#   loop and compares them against the real A4_Lag1PainOnly result, for
#   Best_AUC and each of the four individual model AUCs:
#     - Individual: each person's real value vs. their own null distribution
#       (mean/SD/CI), with an empirical one-sided p-value (proportion of that
#       person's own null draws >= their real value).
#     - Group: consolidate each person's null distribution to their own null
#       mean first, then compare the group's real average against the group
#       average of those null means via a paired t-test.
#
###############################################################################

cat("\n--- A4 null model comparison ---\n")

a4_null_patient_dir <- file.path(results_dir, "A4_NullModel", "PatientResults")
sm_a4_real           <- all_summaries[["A4_Lag1PainOnly"]]

null_metric_cols <- c("Best_AUC", "RandomForest_AUC", "LogisticRegression_AUC",
                      "GaussianProcess_AUC", "Ensemble_AUC")

if (!is.null(sm_a4_real) && "Best_AUC" %in% names(sm_a4_real)) {
  sm_a4_real$Best_AUC <- apply(sm_a4_real[, get_auc_cols(sm_a4_real)], 1, max, na.rm = TRUE)
}

if (is.null(sm_a4_real) || !dir.exists(a4_null_patient_dir)) {
  cat("  A4_Lag1PainOnly summary or A4_NullModel results not found — skipping\n")
} else {

  # Load each participant's null distribution
  null_list <- list()
  for (pid in sm_a4_real$StudyID) {
    null_file <- file.path(a4_null_patient_dir, pid, paste0(pid, "_null_distribution.csv"))
    if (file.exists(null_file))
      null_list[[pid]] <- read_csv(null_file, show_col_types = FALSE) %>%
        select(any_of(c("Repeat", null_metric_cols))) %>%
        mutate(StudyID = pid)
  }

  if (length(null_list) == 0) {
    cat("  No A4_NullModel distributions found on disk — skipping\n")
  } else {
    all_null <- bind_rows(null_list)
    cat("  Null distributions loaded for", length(unique(all_null$StudyID)), "participants\n")

    empirical_p <- function(observed, null_vals)
      (sum(null_vals >= observed) + 1) / (length(null_vals) + 1)

    # Individual comparison, one row per participant x metric
    within_list <- list()
    for (metric in null_metric_cols) {
      if (!(metric %in% names(all_null)) || !(metric %in% names(sm_a4_real))) next
      for (pid in unique(all_null$StudyID)) {
        real_row <- sm_a4_real %>% filter(StudyID == pid)
        if (nrow(real_row) == 0) next
        real_val  <- real_row[[metric]][1]
        null_vals <- all_null[[metric]][all_null$StudyID == pid]
        null_vals <- null_vals[!is.na(null_vals)]
        if (length(null_vals) == 0 || is.na(real_val)) next
        ci    <- calc_ci(null_vals)
        p_raw <- empirical_p(real_val, null_vals)
        within_list[[length(within_list) + 1]] <- data.frame(
          StudyID     = pid,
          Metric      = metric,
          Real_Value  = round(real_val, 3),
          N_Null      = length(null_vals),
          Null_Mean   = round(mean(null_vals), 3),
          Null_SD     = round(sd(null_vals), 3),
          Null_CI_Lo  = round(ci[1], 3),
          Null_CI_Hi  = round(ci[2], 3),
          P_value_raw = p_raw,
          P_value     = ifelse(p_raw < .001, "<0.001", format_decimal(p_raw, 3))
        )
      }
    }
    within_participant <- bind_rows(within_list)

    if (nrow(within_participant) > 0) {
      write_csv(within_participant,
                file.path(analysis_output_dir, "A4_null_within_participant_comparison.csv"))
      cat("  Saved: A4_null_within_participant_comparison.csv\n")
      for (metric in unique(within_participant$Metric)) {
        sub <- within_participant %>% filter(Metric == metric)
        cat("   ", metric, "— participants with p <.05:",
            sum(sub$P_value_raw < .05, na.rm = TRUE), "/", nrow(sub), "\n")
      }
    }

    # Group comparison: consolidate by person first (their own null
    # mean), then paired t-test of Real vs. that consolidated null mean
    # across participants.
    group_list <- list()
    for (metric in unique(within_participant$Metric)) {
      sub <- within_participant %>% filter(Metric == metric)
      if (nrow(sub) < 2) next
      tt  <- t.test(sub$Real_Value, sub$Null_Mean, paired = TRUE)
      dif <- mean(sub$Real_Value - sub$Null_Mean)
      n   <- nrow(sub); se <- sd(sub$Real_Value - sub$Null_Mean) / sqrt(n)
      tc  <- qt(.975, n - 1)
      group_list[[metric]] <- data.frame(
        Metric          = metric,
        N               = n,
        Real_Group_Mean = round(mean(sub$Real_Value), 3),
        Null_Group_Mean = round(mean(sub$Null_Mean), 3),
        Mean_Difference = round(dif, 3),
        CI_Lower        = round(dif - tc*se, 3),
        CI_Upper        = round(dif + tc*se, 3),
        t_statistic     = round(tt$statistic, 2),
        P_value_raw     = tt$p.value,
        P_value         = ifelse(tt$p.value < .001, "<0.001", format_decimal(tt$p.value, 3))
      )
    }
    group_comparison <- bind_rows(group_list)

    if (nrow(group_comparison) > 0) {
      write_csv(group_comparison,
                file.path(analysis_output_dir, "A4_null_group_average_comparison.csv"))
      cat("  Saved: A4_null_group_average_comparison.csv\n")
      print(group_comparison)
    }
  }
}
