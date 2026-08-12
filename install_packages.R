# Target runtime: R 4.5.1
#
# Bundled source dependency:
# mice.reuse.R is included in this repository and sourced by
# step2_imputation_and_modeling.R.
#
# Citation:
# Rockenschaub P. mice.reuse. 2020. Available:
# https://github.com/prockenschaub/Misc/blob/master/R/mice.reuse/mice.reuse.R
# Accessed 2026-08-04.

package_versions <- c(
  "R.utils" = "2.13.0",
  "caret" = "7.0.1",
  "data.table" = "1.17.8",
  "dplyr" = "1.1.4",
  "forcats" = "1.0.0",
  "ggplot2" = "4.0.3",
  "ggridges" = "0.5.6",
  "glmnet" = "4.1-10",
  "kernlab" = "0.9-33",
  "lubridate" = "1.9.4",
  "mice" = "3.18.0",
  "patchwork" = "1.3.2",
  "pROC" = "1.18.5",
  "purrr" = "1.2.1",
  "randomForest" = "4.7-1.2",
  "readr" = "2.1.5",
  "stringr" = "1.5.1",
  "tableone" = "0.13.2",
  "testthat" = "3.2.3",
  "tibble" = "3.3.0",
  "tidyr" = "1.3.1"
)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

for (pkg in names(package_versions)) {
  target_version <- unname(package_versions[pkg])
  is_installed <- requireNamespace(pkg, quietly = TRUE)

  if (nzchar(target_version)) {
    installed_version <- if (is_installed) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }

    if (is.na(installed_version) || installed_version != target_version) {
      message(sprintf("Installing %s %s", pkg, target_version))
      remotes::install_version(pkg, version = target_version, upgrade = "never")
    } else {
      message(sprintf("%s %s already installed", pkg, installed_version))
    }
  } else if (!is_installed) {
    message(sprintf("Installing %s (latest available; fill in version manually if needed)", pkg))
    install.packages(pkg)
  } else {
    message(sprintf("%s already installed", pkg))
  }
}
