# =========================================================
# Script title: 09_models_exp4_field_spm_v2.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Fit and evaluate Experiment 4 field-derived SPM
#          gradient models for motile fraction using primary
#          cell counts, with explicit comparison of simple
#          monotonic versus nonlinear NTU-response forms to
#          test H4: whether NTU alone is an incomplete
#          predictor under environmentally realistic whole-
#          mixture exposures.
# Inputs:
# - data_derived/ms4_analysis_derived.csv or .rds
# Outputs:
# - outputs/models/exp4_field_spm/*.rds
# - outputs/models/exp4_field_spm/*_fixed_effects.csv
# - outputs/models/exp4_field_spm/*_metadata.csv
# - outputs/models/exp4_field_spm/*_diagnostics_summary.csv
# - outputs/models/exp4_field_spm/*_model_comparison.csv
# - outputs/models/exp4_field_spm/*_predictions.csv
# - outputs/tables/09_exp4_field_spm_dataset_summary.csv
# - outputs/tables/09_exp4_field_spm_model_comparison.csv
# - outputs/tables/09_exp4_field_spm_model_summary_table.csv
# - outputs/figures/models_exp4/Fig5_field_spm_model_predictions_*.{pdf,png,tiff}
# - outputs/logs/09_models_exp4_field_spm_log_*.txt
# Date created: 26 March 2026
# Last updated: 30 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 04_derive_variables.R before this script.
# - Primary fitted response:
#   cbind(mobile_cell_count, stationary_cell_count)
# - Manuscript 4 block must be restricted to:
#   experiment_num %in% c(8.2, 9.2, 10.2, 11.2)
# - Focal modelling subset for this script:
#   experiment_num == 8.2
# - Experiment 4 design:
#   Environmentally collected SPM across an NTU gradient,
#   with natural mixed composition.
# - Primary exposure metric for Experiment 4 is NTU.
# - Primary brief model:
#   cbind(mobile_cell_count, stationary_cell_count) ~
#   log10(ntu + 1) * days_from_start + (1 | culture)
# - H4 requires explicit testing of whether a simple
#   monotonic NTU-response is adequate; therefore a
#   nonlinear quadratic alternative is treated as a core
#   comparator rather than an optional side sensitivity.
# - If only one observed day is present, days_from_start is
#   retained as metadata/plot label only and is not fitted.
# - If only one culture level is present, no random effect
#   is included for culture.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 09: MODELS EXP4 FIELD SPM\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check setup objects
# ---------------------------------------------------------

required_objects <- c(
  "project_root",
  "dir_data_derived",
  "dir_tables",
  "dir_logs",
  "project_title",
  "manuscript_short"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "The following required setup object(s) are missing:\n- ",
      paste(missing_objects, collapse = "\n- "),
      "\nPlease run 01_setup_packages_and_paths.R first."
    ),
    call. = FALSE
  )
}

cat("Setup objects verified.\n\n")

# ---------------------------------------------------------
# 2. Source helper scripts if available
# ---------------------------------------------------------

helper_files <- c(
  "helpers_theme.R",
  "helpers_save_figures.R",
  "helpers_labels.R",
  "helpers_tables.R",
  "helpers_model_checks.R"
)

for (hf in helper_files) {
  helper_path <- file.path(project_root, "R", hf)
  if (file.exists(helper_path)) {
    source(helper_path)
    cat("Loaded helper:", hf, "\n")
  }
}
cat("\n")

# ---------------------------------------------------------
# 3. Small local helper functions
# ---------------------------------------------------------

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

write_caption_md_local <- function(figure_name, caption_text, subdir) {
  fig_dir <- file.path(project_root, "outputs", "figures", subdir)
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  caption_file <- file.path(fig_dir, paste0(figure_name, "_caption.md"))
  
  writeLines(
    text = c(
      paste0("# ", figure_name),
      "",
      caption_text
    ),
    con = caption_file
  )
  
  invisible(caption_file)
}

save_rds_local <- function(object, file_path) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = file_path)
  invisible(file_path)
}

rhs_with_random <- function(fixed_rhs, include_random_culture) {
  if (isTRUE(include_random_culture)) {
    paste0(fixed_rhs, " + (1 | culture)")
  } else {
    fixed_rhs
  }
}

as_binomial_formula <- function(rhs) {
  stats::as.formula(
    paste0("cbind(mobile_cell_count, stationary_cell_count) ~ ", rhs)
  )
}

# ---------------------------------------------------------
# 4. Define input and output paths
# ---------------------------------------------------------

file_input_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")

exp4_model_subdir <- "exp4_field_spm"
exp4_figure_subdir <- "models_exp4"

dir_model_output <- file.path(project_root, "outputs", "models", exp4_model_subdir)

file_dataset_summary <- file.path(dir_tables, "09_exp4_field_spm_dataset_summary.csv")
file_model_comparison_table <- file.path(dir_tables, "09_exp4_field_spm_model_comparison.csv")
file_model_summary_table <- file.path(dir_tables, "09_exp4_field_spm_model_summary_table.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_model_log <- file.path(
  dir_logs,
  paste0("09_models_exp4_field_spm_log_", timestamp_now, ".txt")
)

dir.create(dir_model_output, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 5. Import derived dataset
# ---------------------------------------------------------

if (file.exists(file_input_rds)) {
  dat <- readRDS(file_input_rds)
  input_source_used <- file_input_rds
} else if (file.exists(file_input_csv)) {
  dat <- readr::read_csv(file_input_csv, show_col_types = FALSE, progress = FALSE)
  input_source_used <- file_input_csv
} else {
  stop(
    paste0(
      "No derived analysis dataset found.\nExpected one of:\n- ",
      file_input_rds,
      "\n- ",
      file_input_csv,
      "\nPlease run 04_derive_variables.R first."
    ),
    call. = FALSE
  )
}

cat("Derived analysis dataset loaded from:\n")
cat(input_source_used, "\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

# ---------------------------------------------------------
# 6. Check required columns
# ---------------------------------------------------------

required_columns <- c(
  "experiment",
  "experiment_num",
  "experiment_type",
  "mobile_cell_count",
  "stationary_cell_count",
  "ntu",
  "culture"
)

missing_required_columns <- required_columns[!required_columns %in% names(dat)]

if (length(missing_required_columns) > 0) {
  stop(
    paste0(
      "The following required Experiment 4 column(s) are missing:\n- ",
      paste(missing_required_columns, collapse = "\n- "),
      "\nPlease review upstream scripts."
    ),
    call. = FALSE
  )
}

cat("Required Experiment 4 columns verified.\n\n")

# ---------------------------------------------------------
# 7. Restrict to Manuscript 4 block, then focal experiment
# ---------------------------------------------------------

ms4_experiment_nums <- c(8.2, 9.2, 10.2, 11.2)

dat <- dat |>
  dplyr::mutate(
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    experiment_chr = stringr::str_squish(as.character(experiment)),
    experiment_type = as.character(experiment_type)
  ) |>
  dplyr::filter(
    !is.na(experiment_num),
    experiment_num %in% ms4_experiment_nums
  )

cat("Restricted to Manuscript 4 experimental block.\n")
cat("Rows in MS4 block:", nrow(dat), "\n\n")

exp4_dat <- dat |>
  dplyr::filter(experiment_num == 8.2)

if (!"days_from_start" %in% names(exp4_dat)) {
  exp4_dat$days_from_start <- NA_integer_
}

if (nrow(exp4_dat) == 0) {
  stop("No Experiment 4 rows remain after experiment filtering.", call. = FALSE)
}

exp4_dat <- exp4_dat |>
  dplyr::mutate(
    culture = as.factor(culture),
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    ntu = suppressWarnings(as.numeric(ntu))
  ) |>
  dplyr::filter(
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    !is.na(culture),
    !is.na(ntu)
  )

if (nrow(exp4_dat) == 0) {
  stop(
    "No Experiment 4 rows remain after removing missing response/predictor values.",
    call. = FALSE
  )
}

cat("Experiment 4 focal subset created.\n")
cat("Rows:", nrow(exp4_dat), "\n")
cat("Cultures:", dplyr::n_distinct(exp4_dat$culture, na.rm = TRUE), "\n")
cat(
  "Days:",
  if (all(is.na(exp4_dat$days_from_start))) {
    "none available"
  } else {
    paste(sort(unique(stats::na.omit(exp4_dat$days_from_start))), collapse = ", ")
  },
  "\n"
)
cat("NTU range:", paste0(round(safe_min(exp4_dat$ntu), 3), " to ", round(safe_max(exp4_dat$ntu), 3)), "\n\n")

# ---------------------------------------------------------
# 8. Derive support variables
# ---------------------------------------------------------

exp4_dat <- exp4_dat |>
  dplyr::mutate(
    total_cells = mobile_cell_count + stationary_cell_count,
    motility_ratio = dplyr::if_else(
      total_cells > 0,
      mobile_cell_count / total_cells,
      NA_real_
    ),
    log10_ntu_plus1 = log10(ntu + 1),
    ntu_sq = log10_ntu_plus1^2,
    day_label = dplyr::if_else(
      !is.na(days_from_start),
      paste0("Day ", days_from_start),
      "Single sampling day"
    ),
    experiment_plot_label = "Experiment 4: Field-derived SPM"
  ) |>
  dplyr::filter(total_cells > 0)

if (nrow(exp4_dat) == 0) {
  stop("No valid Experiment 4 rows remain after deriving support variables.", call. = FALSE)
}

if (sum(!is.na(exp4_dat$log10_ntu_plus1)) < 2) {
  stop(
    "Experiment 4 modelling requires non-missing NTU values.",
    call. = FALSE
  )
}

observed_days <- sort(unique(stats::na.omit(exp4_dat$days_from_start)))
n_day_levels <- length(observed_days)
n_culture_levels <- dplyr::n_distinct(exp4_dat$culture, na.rm = TRUE)

has_day_variation <- n_day_levels > 1
has_culture_variation <- n_culture_levels > 1

cat("Experiment 4 modelling support variables derived.\n")
cat("Day variation available:", has_day_variation, "\n")
cat("Culture variation available:", has_culture_variation, "\n\n")

# ---------------------------------------------------------
# 9. Build dataset summary table
# ---------------------------------------------------------

has_well <- "well" %in% names(exp4_dat)
has_video_file <- "video_file" %in% names(exp4_dat)

group_vars <- character(0)
if (!all(is.na(exp4_dat$days_from_start))) {
  group_vars <- "days_from_start"
}

if (length(group_vars) > 0) {
  exp4_dataset_summary <- exp4_dat |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
      n_wells = if (has_well) dplyr::n_distinct(well, na.rm = TRUE) else NA_integer_,
      n_videos = if (has_video_file) dplyr::n_distinct(video_file, na.rm = TRUE) else NA_integer_,
      min_ntu = safe_min(ntu),
      max_ntu = safe_max(ntu),
      mean_ntu = safe_mean(ntu),
      mean_total_cells = safe_mean(total_cells),
      mean_motility_ratio = safe_mean(motility_ratio),
      sd_motility_ratio = safe_sd(motility_ratio),
      .groups = "drop"
    )
} else {
  exp4_dataset_summary <- exp4_dat |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
      n_wells = if (has_well) dplyr::n_distinct(well, na.rm = TRUE) else NA_integer_,
      n_videos = if (has_video_file) dplyr::n_distinct(video_file, na.rm = TRUE) else NA_integer_,
      min_ntu = safe_min(ntu),
      max_ntu = safe_max(ntu),
      mean_ntu = safe_mean(ntu),
      mean_total_cells = safe_mean(total_cells),
      mean_motility_ratio = safe_mean(motility_ratio),
      sd_motility_ratio = safe_sd(motility_ratio)
    )
}

readr::write_csv(exp4_dataset_summary, file_dataset_summary)

cat("Experiment 4 dataset summary written to:\n")
cat(file_dataset_summary, "\n\n")

# ---------------------------------------------------------
# 10. Define candidate model formulas
# ---------------------------------------------------------

candidate_formula_strings <- list()

if (has_day_variation) {
  candidate_formula_strings$monotonic_primary <- rhs_with_random(
    "log10_ntu_plus1 * days_from_start",
    has_culture_variation
  )
  
  candidate_formula_strings$monotonic_additive <- rhs_with_random(
    "log10_ntu_plus1 + days_from_start",
    has_culture_variation
  )
  
  candidate_formula_strings$nonlinear_primary <- rhs_with_random(
    "log10_ntu_plus1 + I(log10_ntu_plus1^2) + days_from_start",
    has_culture_variation
  )
  
  candidate_formula_strings$day_only <- rhs_with_random(
    "days_from_start",
    has_culture_variation
  )
} else {
  candidate_formula_strings$monotonic_primary <- rhs_with_random(
    "log10_ntu_plus1",
    has_culture_variation
  )
  
  candidate_formula_strings$nonlinear_primary <- rhs_with_random(
    "log10_ntu_plus1 + I(log10_ntu_plus1^2)",
    has_culture_variation
  )
}

candidate_formula_strings$null <- rhs_with_random(
  "1",
  has_culture_variation
)

candidate_formulas <- lapply(candidate_formula_strings, as_binomial_formula)

cat("Candidate model formulas defined.\n")
for (nm in names(candidate_formula_strings)) {
  cat("-", nm, ":", candidate_formula_strings[[nm]], "\n")
}
cat("\n")

# ---------------------------------------------------------
# 11. Fit candidate models
# ---------------------------------------------------------

if (!exists("fit_check_save_model", mode = "function", inherits = TRUE)) {
  stop("helpers_model_checks.R was not loaded correctly.", call. = FALSE)
}

fit_results <- list()
failed_models <- character(0)

for (nm in names(candidate_formulas)) {
  fit_results[[nm]] <- tryCatch(
    fit_check_save_model(
      formula = candidate_formulas[[nm]],
      data = exp4_dat,
      model_name = paste0("exp4_field_spm_", nm),
      model_subdir = exp4_model_subdir,
      save_terms_csv = TRUE,
      save_meta_csv = TRUE,
      save_diagnostics = TRUE,
      run_dharma = TRUE,
      run_performance = TRUE,
      conf.level = 0.95,
      exponentiate = FALSE,
      quiet = FALSE
    ),
    error = function(e) {
      message("[model] Model failed: ", nm, " | ", conditionMessage(e))
      failed_models <<- c(failed_models, nm)
      NULL
    }
  )
}

fit_results <- fit_results[!vapply(fit_results, is.null, logical(1))]

cat("Candidate model fitting complete.\n")
if (length(failed_models) > 0) {
  cat("Failed candidate models:", paste(failed_models, collapse = ", "), "\n")
}
cat("\n")

if (length(fit_results) < 2) {
  stop(
    paste0(
      "Fewer than two candidate models fitted successfully.\n",
      "Failed models: ",
      ifelse(length(failed_models) == 0, "none recorded", paste(failed_models, collapse = ", "))
    ),
    call. = FALSE
  )
}

cat("At least two candidate models fitted successfully.\n\n")

# ---------------------------------------------------------
# 12. Compare candidate models
# ---------------------------------------------------------

if (!exists("compare_models_aic", mode = "function", inherits = TRUE)) {
  stop("compare_models_aic() is not available.", call. = FALSE)
}

comparison_input <- lapply(fit_results, function(x) x$model)

model_comparison <- do.call(compare_models_aic, comparison_input) |>
  dplyr::mutate(
    experiment = "Experiment 4",
    experiment_num = 8.2,
    complexity_rank = dplyr::case_when(
      model == "nonlinear_primary" ~ 4L,
      model == "monotonic_primary" ~ 3L,
      model == "monotonic_additive" ~ 2L,
      model == "day_only" ~ 1L,
      model == "null" ~ 0L,
      TRUE ~ 99L
    )
  ) |>
  dplyr::arrange(delta_aic, dplyr::desc(complexity_rank))

readr::write_csv(model_comparison, file_model_comparison_table)

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = model_comparison,
    file_stem = "exp4_field_spm_model_comparison",
    subdir = exp4_model_subdir,
    quiet = FALSE
  )
}

cat("Model comparison written to:\n")
cat(file_model_comparison_table, "\n\n")

# ---------------------------------------------------------
# 13. Likelihood ratio comparisons
# ---------------------------------------------------------

lrt_monotonic_additive_vs_primary <- NULL
lrt_day_only_vs_monotonic_additive <- NULL
lrt_null_vs_day_only <- NULL

if (all(c("monotonic_additive", "monotonic_primary") %in% names(fit_results))) {
  lrt_monotonic_additive_vs_primary <- tryCatch(
    stats::anova(fit_results$monotonic_additive$model, fit_results$monotonic_primary$model),
    error = function(e) e
  )
}

if (all(c("day_only", "monotonic_additive") %in% names(fit_results))) {
  lrt_day_only_vs_monotonic_additive <- tryCatch(
    stats::anova(fit_results$day_only$model, fit_results$monotonic_additive$model),
    error = function(e) e
  )
}

if (all(c("null", "day_only") %in% names(fit_results))) {
  lrt_null_vs_day_only <- tryCatch(
    stats::anova(fit_results$null$model, fit_results$day_only$model),
    error = function(e) e
  )
}

# Direct H4 support note
has_monotonic_model <- "monotonic_primary" %in% model_comparison$model
has_nonlinear_model <- "nonlinear_primary" %in% model_comparison$model

h4_support_note <- if (has_monotonic_model && has_nonlinear_model) {
  mono_aic <- model_comparison$AIC[model_comparison$model == "monotonic_primary"][1]
  nonlin_aic <- model_comparison$AIC[model_comparison$model == "nonlinear_primary"][1]
  delta_nonlin_vs_mono <- mono_aic - nonlin_aic
  
  if (is.na(delta_nonlin_vs_mono)) {
    "H4 comparison note unavailable due to missing AIC values."
  } else if (delta_nonlin_vs_mono > 2) {
    paste0(
      "The nonlinear NTU model outperformed the monotonic NTU model by ",
      round(delta_nonlin_vs_mono, 3),
      " AIC units, supporting the interpretation that a simple monotonic NTU-response is inadequate for field-derived SPM."
    )
  } else if (delta_nonlin_vs_mono >= 0 && delta_nonlin_vs_mono <= 2) {
    paste0(
      "The nonlinear NTU model was within ",
      round(delta_nonlin_vs_mono, 3),
      " AIC units of the monotonic NTU model, suggesting that simple monotonic NTU-response may be insufficient and should be interpreted cautiously."
    )
  } else {
    paste0(
      "The monotonic NTU model outperformed the nonlinear NTU model by ",
      round(abs(delta_nonlin_vs_mono), 3),
      " AIC units."
    )
  }
} else {
  "Both monotonic and nonlinear core Experiment 4 models were not available for direct H4 comparison."
}

# ---------------------------------------------------------
# 14. Select preferred model
# ---------------------------------------------------------

best_delta <- min(model_comparison$delta_aic, na.rm = TRUE)

preferred_model_name <- model_comparison |>
  dplyr::filter(delta_aic == best_delta) |>
  dplyr::arrange(dplyr::desc(complexity_rank)) |>
  dplyr::slice(1) |>
  dplyr::pull(model)

preferred_model <- fit_results[[preferred_model_name]]$model

cat("Preferred model selected:\n")
cat(preferred_model_name, "\n\n")

# ---------------------------------------------------------
# 15. Build prediction grid
# ---------------------------------------------------------

if (!exists("make_prediction_grid", mode = "function", inherits = TRUE) ||
    !exists("predict_glmmtmb_response", mode = "function", inherits = TRUE)) {
  stop("Prediction helpers are not available.", call. = FALSE)
}

pred_at <- list()
if (has_day_variation) {
  pred_at$days_from_start <- sort(unique(exp4_dat$days_from_start))
}

if (preferred_model_name %in% c("monotonic_primary", "monotonic_additive", "nonlinear_primary")) {
  pred_grid <- make_prediction_grid(
    data = exp4_dat,
    focal_terms = c("log10_ntu_plus1"),
    at = pred_at,
    n_numeric = 100
  )
  
  pred_out <- predict_glmmtmb_response(
    model = preferred_model,
    newdata = pred_grid,
    conf.level = 0.95,
    re.form = NA
  ) |>
    dplyr::mutate(
      ntu_backtransformed = (10^log10_ntu_plus1) - 1,
      day_label = if ("days_from_start" %in% names(pred_grid)) {
        dplyr::if_else(!is.na(days_from_start), paste0("Day ", days_from_start), "Day 4")
      } else {
        "Day 4"
      },
      prediction_scale = "log10_ntu_plus1"
    )
} else if (preferred_model_name == "day_only") {
  pred_grid <- make_prediction_grid(
    data = exp4_dat,
    focal_terms = c("days_from_start"),
    at = list(
      days_from_start = sort(unique(exp4_dat$days_from_start)),
      log10_ntu_plus1 = stats::median(exp4_dat$log10_ntu_plus1, na.rm = TRUE)
    ),
    n_numeric = 100
  )
  
  pred_out <- predict_glmmtmb_response(
    model = preferred_model,
    newdata = pred_grid,
    conf.level = 0.95,
    re.form = NA
  ) |>
    dplyr::mutate(
      ntu_backtransformed = (10^log10_ntu_plus1) - 1,
      day_label = paste0("Day ", days_from_start),
      prediction_scale = "day_only"
    )
} else {
  pred_grid <- make_prediction_grid(
    data = exp4_dat,
    focal_terms = c("days_from_start"),
    at = list(
      days_from_start = if (has_day_variation) sort(unique(exp4_dat$days_from_start)) else NA_integer_,
      log10_ntu_plus1 = stats::median(exp4_dat$log10_ntu_plus1, na.rm = TRUE)
    ),
    n_numeric = 100
  )
  
  pred_out <- predict_glmmtmb_response(
    model = preferred_model,
    newdata = pred_grid,
    conf.level = 0.95,
    re.form = NA
  ) |>
    dplyr::mutate(
      ntu_backtransformed = (10^log10_ntu_plus1) - 1,
      day_label = if ("days_from_start" %in% names(pred_grid)) paste0("Day ", days_from_start) else "Day 4",
      prediction_scale = "null_model"
    )
}

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = pred_out,
    file_stem = "exp4_field_spm_preferred_model_predictions",
    subdir = exp4_model_subdir,
    quiet = FALSE
  )
}

cat("Prediction outputs generated and saved.\n\n")

# ---------------------------------------------------------
# 16. Build fixed-effect summary for preferred model
# ---------------------------------------------------------

if (!exists("build_model_term_summary", mode = "function", inherits = TRUE)) {
  stop("build_model_term_summary() is not available.", call. = FALSE)
}

preferred_fixed_effects <- tryCatch(
  build_model_term_summary(
    model = preferred_model,
    model_name = preferred_model_name,
    conf.level = 0.95,
    exponentiate = FALSE
  ) |>
    dplyr::mutate(
      experiment = "Experiment 4",
      experiment_num = 8.2
    ),
  error = function(e) {
    tibble::tibble(
      model_name = preferred_model_name,
      term = NA_character_,
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      note = paste("Fixed-effect summary not available for selected model:", e$message),
      experiment = "Experiment 4",
      experiment_num = 8.2
    )
  }
)

readr::write_csv(preferred_fixed_effects, file_model_summary_table)

cat("Preferred model fixed-effect summary written to:\n")
cat(file_model_summary_table, "\n\n")

# ---------------------------------------------------------
# 17. Build Experiment 4 figure
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

if (exists("particle_palette", mode = "function", inherits = TRUE)) {
  exp4_particle_palette <- particle_palette()
} else if (exists("particle_palette", inherits = TRUE)) {
  exp4_particle_palette <- get("particle_palette", inherits = TRUE)
} else {
  exp4_particle_palette <- c("SPM" = "#0072B2")
}

exp4_spm_colour <- unname(if ("SPM" %in% names(exp4_particle_palette)) exp4_particle_palette["SPM"] else "#0072B2")

plot_raw <- exp4_dat |>
  dplyr::mutate(
    day_label = dplyr::if_else(
      !is.na(days_from_start),
      paste0("Day ", days_from_start),
      "Day 4"
    )
  )

if (preferred_model_name %in% c("monotonic_primary", "monotonic_additive", "nonlinear_primary")) {
  p_exp4 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_raw,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = motility_ratio
      ),
      alpha = 0.30,
      size = 1.8,
      colour = exp4_spm_colour
    ) +
    ggplot2::geom_ribbon(
      data = pred_out,
      ggplot2::aes(
        x = log10_ntu_plus1,
        ymin = conf_low_response,
        ymax = conf_high_response
      ),
      alpha = 0.18,
      fill = exp4_spm_colour
    ) +
    ggplot2::geom_line(
      data = pred_out,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = fit_response
      ),
      linewidth = 0.8,
      colour = exp4_spm_colour
    ) +
    ggplot2::labs(
      title = "Field-derived SPM gradient",
      subtitle = if (preferred_model_name == "nonlinear_primary") {
        "Model-predicted motile fraction versus log10(NTU + 1): nonlinear NTU-response"
      } else {
        "Model-predicted motile fraction versus log10(NTU + 1)"
      },
      x = "log10(NTU + 1)",
      y = "Motile fraction",
      caption = if (has_day_variation) {
        "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals. Panels are shown by day from start."
      } else {
        "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals. Experiment 4 was sampled on a single day only, so no day facet is shown."
      }
    )
  
  if (has_day_variation) {
    p_exp4 <- p_exp4 + ggplot2::facet_wrap(~ day_label)
  }
} else {
  p_exp4 <- ggplot2::ggplot() +
    ggplot2::geom_boxplot(
      data = plot_raw,
      ggplot2::aes(
        x = factor(day_label),
        y = motility_ratio
      ),
      outlier.alpha = 0.35,
      width = 0.6,
      colour = exp4_spm_colour
    ) +
    ggplot2::geom_jitter(
      data = plot_raw,
      ggplot2::aes(
        x = factor(day_label),
        y = motility_ratio
      ),
      width = 0.12,
      alpha = 0.25,
      size = 1.4,
      colour = exp4_spm_colour
    ) +
    ggplot2::geom_pointrange(
      data = pred_out,
      ggplot2::aes(
        x = factor(day_label),
        y = fit_response,
        ymin = conf_low_response,
        ymax = conf_high_response
      ),
      linewidth = 0.6,
      colour = "black"
    ) +
    ggplot2::labs(
      title = "Field-derived SPM gradient",
      subtitle = "Day-based fitted values from preferred model",
      x = "Sampling day",
      y = "Motile fraction",
      caption = "Points show raw observations and pointranges show model predictions with 95% confidence intervals."
    )
}

# ---------------------------------------------------------
# 18. Save Experiment 4 figure
# ---------------------------------------------------------

saved_figures <- character(0)

if (exists("save_figure_both_widths", mode = "function", inherits = TRUE)) {
  saved_figures <- save_figure_both_widths(
    plot = p_exp4,
    figure_name = "Fig5_field_spm_model_predictions",
    subdir = exp4_figure_subdir,
    height = "standard",
    quiet = TRUE
  )
}

caption_text <- if (preferred_model_name == "nonlinear_primary") {
  "Model-predicted motile fraction as a function of log10(NTU + 1) in the field-derived SPM gradient experiment. A nonlinear NTU-response model was selected to evaluate whether simple monotonic NTU alone adequately predicts motility. Raw observations are overlaid, 95% confidence intervals are shown, and panels are facetted by day from start where appropriate."
} else {
  "Model-predicted motile fraction as a function of log10(NTU + 1) in the field-derived SPM gradient experiment. Raw observations are overlaid, 95% confidence intervals are shown, and panels are facetted by day from start where appropriate."
}

if (exists("write_figure_caption_md", mode = "function", inherits = TRUE)) {
  write_figure_caption_md(
    figure_name = "Fig5_field_spm_model_predictions",
    caption_text = caption_text,
    subdir = exp4_figure_subdir
  )
} else {
  write_caption_md_local(
    figure_name = "Fig5_field_spm_model_predictions",
    caption_text = caption_text,
    subdir = exp4_figure_subdir
  )
}

cat("Experiment 4 figure processed.\n\n")

# ---------------------------------------------------------
# 19. Save extra model objects and summaries
# ---------------------------------------------------------

selection_bundle <- list(
  observed_days = observed_days,
  n_day_levels = n_day_levels,
  n_culture_levels = n_culture_levels,
  has_day_variation = has_day_variation,
  has_culture_variation = has_culture_variation,
  candidate_formula_strings = candidate_formula_strings,
  failed_models = failed_models,
  lrt_monotonic_additive_vs_primary = lrt_monotonic_additive_vs_primary,
  lrt_day_only_vs_monotonic_additive = lrt_day_only_vs_monotonic_additive,
  lrt_null_vs_day_only = lrt_null_vs_day_only,
  h4_support_note = h4_support_note,
  preferred_model_name = preferred_model_name
)

if (exists("write_model_summary_rds", mode = "function", inherits = TRUE)) {
  write_model_summary_rds(
    object = selection_bundle,
    file_stem = "exp4_field_spm_model_tests_and_selection",
    subdir = exp4_model_subdir,
    quiet = FALSE
  )
  
  write_model_summary_rds(
    object = pred_out,
    file_stem = "exp4_field_spm_preferred_model_predictions",
    subdir = exp4_model_subdir,
    quiet = FALSE
  )
} else {
  save_rds_local(
    object = selection_bundle,
    file_path = file.path(dir_model_output, "exp4_field_spm_model_tests_and_selection.rds")
  )
  
  save_rds_local(
    object = pred_out,
    file_path = file.path(dir_model_output, "exp4_field_spm_preferred_model_predictions.rds")
  )
}

# ---------------------------------------------------------
# 20. Optional emmeans summaries
# ---------------------------------------------------------

emmeans_day <- NULL

if (exists("get_emmeans_table", mode = "function", inherits = TRUE) &&
    has_day_variation &&
    preferred_model_name %in% c("day_only", "monotonic_additive")) {
  
  emmeans_day <- tryCatch(
    get_emmeans_table(
      model = preferred_model,
      specs = ~ days_from_start,
      type = "response"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(emmeans_day) && exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
    write_model_summary_csv(
      data = emmeans_day,
      file_stem = "exp4_field_spm_emmeans_by_day",
      subdir = exp4_model_subdir,
      quiet = FALSE
    )
  }
}

# ---------------------------------------------------------
# 21. Write log
# ---------------------------------------------------------

sink(file_model_log)
cat("SPM Analysis - 09_models_exp4_field_spm log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")

cat("Input dataset:\n")
cat(input_source_used, "\n\n")

cat("Experiment 4 subset dimensions:\n")
cat("Rows:", nrow(exp4_dat), "\n")
cat("Columns:", ncol(exp4_dat), "\n\n")

cat("Observed days:\n")
if (length(observed_days) == 0) {
  cat("No non-missing day values available.\n\n")
} else {
  cat(paste(observed_days, collapse = ", "), "\n\n")
}

cat("Observed NTU range:\n")
cat(paste0(round(safe_min(exp4_dat$ntu), 3), " to ", round(safe_max(exp4_dat$ntu), 3)), "\n\n")

cat("Detected design features:\n")
cat("has_day_variation:", has_day_variation, "\n")
cat("has_culture_variation:", has_culture_variation, "\n")
cat("n_culture_levels:", n_culture_levels, "\n\n")

cat("Candidate formulas used:\n")
for (nm in names(candidate_formula_strings)) {
  cat("-", nm, ":", candidate_formula_strings[[nm]], "\n")
}
cat("\n")

if (length(failed_models) > 0) {
  cat("Failed candidate models:\n")
  cat(paste(failed_models, collapse = ", "), "\n\n")
}

cat("Experiment 4 dataset summary:\n")
print(exp4_dataset_summary)
cat("\n")

cat("Model comparison:\n")
print(model_comparison)
cat("\n")

cat("Preferred model name:\n")
cat(preferred_model_name, "\n\n")

cat("Preferred fixed effects:\n")
print(preferred_fixed_effects)
cat("\n")

cat("Likelihood ratio test: monotonic additive vs monotonic primary\n")
print(lrt_monotonic_additive_vs_primary)
cat("\n")

cat("Likelihood ratio test: day-only vs monotonic additive\n")
print(lrt_day_only_vs_monotonic_additive)
cat("\n")

cat("Likelihood ratio test: null vs day-only\n")
print(lrt_null_vs_day_only)
cat("\n")

cat("H4 support note:\n")
cat(h4_support_note, "\n\n")

if (!is.null(emmeans_day)) {
  cat("Estimated marginal means by day:\n")
  print(emmeans_day)
  cat("\n")
}

cat("Saved figures:\n")
if (length(saved_figures) == 0) {
  cat("No figures were automatically saved.\n")
} else {
  cat(paste(saved_figures, collapse = "\n"), "\n")
}
cat("\n")

cat("Session information:\n\n")
print(utils::sessionInfo())
sink()

cat("Model log written to:\n")
cat(file_model_log, "\n\n")

# ---------------------------------------------------------
# 23. Sensitivity analysis: Day 4 only
# ---------------------------------------------------------

cat("Running Experiment 4 Day 4 sensitivity analysis...\n")

exp4_day4_dat <- exp4_dat |>
  dplyr::filter(days_from_start == 4)

if (nrow(exp4_day4_dat) == 0) {
  warning("No Day 4 rows available for Experiment 4 sensitivity analysis.", call. = FALSE)
} else {
  
  cat("Day 4 sensitivity subset created.\n")
  cat("Rows:", nrow(exp4_day4_dat), "\n")
  cat("Cultures:", dplyr::n_distinct(exp4_day4_dat$culture, na.rm = TRUE), "\n")
  cat(
    "NTU range:",
    paste0(
      round(safe_min(exp4_day4_dat$ntu), 3),
      " to ",
      round(safe_max(exp4_day4_dat$ntu), 3)
    ),
    "\n\n"
  )
  
  has_day4_culture_variation <- dplyr::n_distinct(exp4_day4_dat$culture, na.rm = TRUE) > 1
  
  day4_formula_strings <- list(
    monotonic_primary = rhs_with_random(
      "log10_ntu_plus1",
      has_day4_culture_variation
    ),
    nonlinear_primary = rhs_with_random(
      "log10_ntu_plus1 + I(log10_ntu_plus1^2)",
      has_day4_culture_variation
    ),
    null = rhs_with_random(
      "1",
      has_day4_culture_variation
    )
  )
  
  day4_formulas <- lapply(day4_formula_strings, as_binomial_formula)
  
  day4_fit_results <- list()
  day4_failed_models <- character(0)
  
  for (nm in names(day4_formulas)) {
    day4_fit_results[[nm]] <- tryCatch(
      fit_check_save_model(
        formula = day4_formulas[[nm]],
        data = exp4_day4_dat,
        model_name = paste0("exp4_field_spm_day4_", nm),
        model_subdir = exp4_model_subdir,
        save_terms_csv = TRUE,
        save_meta_csv = TRUE,
        save_diagnostics = TRUE,
        run_dharma = TRUE,
        run_performance = TRUE,
        conf.level = 0.95,
        exponentiate = FALSE,
        quiet = FALSE
      ),
      error = function(e) {
        message("[model] Day 4 sensitivity model failed: ", nm, " | ", conditionMessage(e))
        day4_failed_models <<- c(day4_failed_models, nm)
        NULL
      }
    )
  }
  
  day4_fit_results <- day4_fit_results[!vapply(day4_fit_results, is.null, logical(1))]
  
  if (length(day4_fit_results) >= 2) {
    
    day4_comparison_input <- lapply(day4_fit_results, function(x) x$model)
    
    day4_model_comparison <- do.call(compare_models_aic, day4_comparison_input) |>
      dplyr::mutate(
        experiment = "Experiment 4 Day 4 sensitivity",
        experiment_num = 8.2,
        subset = "Day 4",
        complexity_rank = dplyr::case_when(
          model == "nonlinear_primary" ~ 3L,
          model == "monotonic_primary" ~ 2L,
          model == "null" ~ 1L,
          TRUE ~ 99L
        )
      ) |>
      dplyr::arrange(delta_aic, dplyr::desc(complexity_rank))
    
    if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
      write_model_summary_csv(
        data = day4_model_comparison,
        file_stem = "exp4_field_spm_day4_model_comparison",
        subdir = exp4_model_subdir,
        quiet = FALSE
      )
    }
    
    day4_best_delta <- min(day4_model_comparison$delta_aic, na.rm = TRUE)
    
    day4_preferred_model_name <- day4_model_comparison |>
      dplyr::filter(delta_aic == day4_best_delta) |>
      dplyr::arrange(dplyr::desc(complexity_rank)) |>
      dplyr::slice(1) |>
      dplyr::pull(model)
    
    day4_preferred_model <- day4_fit_results[[day4_preferred_model_name]]$model
    
    day4_pred_grid <- make_prediction_grid(
      data = exp4_day4_dat,
      focal_terms = c("log10_ntu_plus1"),
      at = list(),
      n_numeric = 100
    )
    
    day4_pred_out <- predict_glmmtmb_response(
      model = day4_preferred_model,
      newdata = day4_pred_grid,
      conf.level = 0.95,
      re.form = NA
    ) |>
      dplyr::mutate(
        ntu_backtransformed = (10^log10_ntu_plus1) - 1,
        day_label = "Day 4",
        prediction_scale = "log10_ntu_plus1",
        subset = "Day 4"
      )
    
    if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
      write_model_summary_csv(
        data = day4_pred_out,
        file_stem = "exp4_field_spm_day4_preferred_model_predictions",
        subdir = exp4_model_subdir,
        quiet = FALSE
      )
    }
    
    day4_fixed_effects <- tryCatch(
      build_model_term_summary(
        model = day4_preferred_model,
        model_name = paste0("day4_", day4_preferred_model_name),
        conf.level = 0.95,
        exponentiate = FALSE
      ) |>
        dplyr::mutate(
          experiment = "Experiment 4 Day 4 sensitivity",
          experiment_num = 8.2,
          subset = "Day 4"
        ),
      error = function(e) {
        tibble::tibble(
          model_name = paste0("day4_", day4_preferred_model_name),
          term = NA_character_,
          estimate = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          note = paste("Fixed-effect summary not available for Day 4 model:", e$message),
          experiment = "Experiment 4 Day 4 sensitivity",
          experiment_num = 8.2,
          subset = "Day 4"
        )
      }
    )
    
    if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
      write_model_summary_csv(
        data = day4_fixed_effects,
        file_stem = "exp4_field_spm_day4_model_summary_table",
        subdir = exp4_model_subdir,
        quiet = FALSE
      )
    }
    
    day4_h4_note <- if (all(c("monotonic_primary", "nonlinear_primary") %in% day4_model_comparison$model)) {
      mono_aic_day4 <- day4_model_comparison$AIC[day4_model_comparison$model == "monotonic_primary"][1]
      nonlin_aic_day4 <- day4_model_comparison$AIC[day4_model_comparison$model == "nonlinear_primary"][1]
      delta_day4 <- mono_aic_day4 - nonlin_aic_day4
      
      if (delta_day4 > 2) {
        paste0(
          "Day 4 sensitivity analysis supported a nonlinear NTU-response over a monotonic NTU-response by ",
          round(delta_day4, 3),
          " AIC units."
        )
      } else if (delta_day4 >= 0) {
        paste0(
          "Day 4 sensitivity analysis gave similar support to nonlinear and monotonic NTU models (ΔAIC = ",
          round(delta_day4, 3),
          ")."
        )
      } else {
        paste0(
          "Day 4 sensitivity analysis favoured the monotonic NTU model by ",
          round(abs(delta_day4), 3),
          " AIC units."
        )
      }
    } else {
      "Day 4 sensitivity analysis did not fit both monotonic and nonlinear core models."
    }
    
    p_exp4_day4 <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = exp4_day4_dat,
        ggplot2::aes(
          x = log10_ntu_plus1,
          y = motility_ratio
        ),
        alpha = 0.30,
        size = 1.8,
        colour = exp4_spm_colour
      ) +
      ggplot2::geom_ribbon(
        data = day4_pred_out,
        ggplot2::aes(
          x = log10_ntu_plus1,
          ymin = conf_low_response,
          ymax = conf_high_response
        ),
        alpha = 0.18,
        fill = exp4_spm_colour
      ) +
      ggplot2::geom_line(
        data = day4_pred_out,
        ggplot2::aes(
          x = log10_ntu_plus1,
          y = fit_response
        ),
        linewidth = 0.8,
        colour = exp4_spm_colour
      ) +
      ggplot2::labs(
        title = "Field-derived SPM gradient: Day 4 sensitivity analysis",
        subtitle = if (day4_preferred_model_name == "nonlinear_primary") {
          "Model-predicted motile fraction versus log10(NTU + 1): nonlinear NTU-response"
        } else {
          "Model-predicted motile fraction versus log10(NTU + 1)"
        },
        x = "log10(NTU + 1)",
        y = "Motile fraction",
        caption = "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals for the Day 4 subset only."
      )
    
    day4_saved_figures <- character(0)
    
    if (exists("save_figure_both_widths", mode = "function", inherits = TRUE)) {
      day4_saved_figures <- save_figure_both_widths(
        plot = p_exp4_day4,
        figure_name = "FigSx_field_spm_day4_sensitivity",
        subdir = exp4_figure_subdir,
        height = "standard",
        quiet = TRUE
      )
    }
    
    day4_caption_text <- "Model-predicted motile fraction as a function of log10(NTU + 1) for the Day 4 subset of the field-derived SPM experiment. Raw observations are overlaid, and 95% confidence intervals are shown. This sensitivity analysis was used to assess whether the nonlinear turbidity–motility relationship identified in the full Experiment 4 dataset was already evident at the earliest sampling point."
    
    if (exists("write_figure_caption_md", mode = "function", inherits = TRUE)) {
      write_figure_caption_md(
        figure_name = "FigSx_field_spm_day4_sensitivity",
        caption_text = day4_caption_text,
        subdir = exp4_figure_subdir
      )
    } else {
      write_caption_md_local(
        figure_name = "FigSx_field_spm_day4_sensitivity",
        caption_text = day4_caption_text,
        subdir = exp4_figure_subdir
      )
    }
    
    day4_bundle <- list(
      n_rows = nrow(exp4_day4_dat),
      n_cultures = dplyr::n_distinct(exp4_day4_dat$culture, na.rm = TRUE),
      ntu_range = c(safe_min(exp4_day4_dat$ntu), safe_max(exp4_day4_dat$ntu)),
      candidate_formula_strings = day4_formula_strings,
      failed_models = day4_failed_models,
      preferred_model_name = day4_preferred_model_name,
      h4_support_note = day4_h4_note
    )
    
    if (exists("write_model_summary_rds", mode = "function", inherits = TRUE)) {
      write_model_summary_rds(
        object = day4_bundle,
        file_stem = "exp4_field_spm_day4_model_tests_and_selection",
        subdir = exp4_model_subdir,
        quiet = FALSE
      )
    } else {
      save_rds_local(
        object = day4_bundle,
        file_path = file.path(dir_model_output, "exp4_field_spm_day4_model_tests_and_selection.rds")
      )
    }
    
    cat("Day 4 sensitivity model comparison:\n")
    print(day4_model_comparison)
    cat("\n")
    
    cat("Day 4 preferred model:\n")
    cat(day4_preferred_model_name, "\n\n")
    
    cat("Day 4 H4 note:\n")
    cat(day4_h4_note, "\n\n")
    
  } else {
    warning(
      "Fewer than two Day 4 sensitivity models fitted successfully; sensitivity comparison not completed.",
      call. = FALSE
    )
  }
}
# ---------------------------------------------------------
# 24. Console summary
# ---------------------------------------------------------

cat("Model comparison:\n")
print(model_comparison)
cat("\n")

cat("Preferred model:\n")
cat(preferred_model_name, "\n\n")

cat("Preferred fixed effects:\n")
print(preferred_fixed_effects)
cat("\n")

cat("H4 support note:\n")
cat(h4_support_note, "\n\n")

cat("09 Experiment 4 modelling outputs written successfully.\n\n")

cat("========================================================\n")
cat("SCRIPT 09 COMPLETE: MODELS EXP4 FIELD SPM\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

