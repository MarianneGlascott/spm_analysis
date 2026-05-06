# =========================================================
# Script title: 06b_models_exp1_light_day4_01.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Manuscript: Manuscript 4 / Article 4
# Purpose: Fit Day 4-only Experiment 1 light-gradient models
#          for kelp zoospore motility using primary count data.
#
# Inputs:
# - data_derived/day4/ms4_day4_analysis_derived.rds or .csv
#
# Outputs:
# - outputs/models/day4/exp1_light_day4/*.rds
# - outputs/models/day4/exp1_light_day4/*.csv
# - outputs/tables/day4/06b_exp1_light_day4_dataset_summary.csv
# - outputs/tables/day4/06b_exp1_light_day4_model_comparison.csv
# - outputs/tables/day4/06b_exp1_light_day4_predictions.csv
# - outputs/tables/day4/06b_exp1_light_day4_relative_effects.csv
# - outputs/figures/day4/Fig4_exp1_light_day4_model_predictions.{pdf,png,tiff}
# - outputs/logs/day4/06b_exp1_light_day4_log_*.txt
#
# Notes:
# - This is the Day 4 analysis branch.
# - This script does not modify the full time-course analysis.
# - Primary response:
#   cbind(mobile_cell_count, stationary_cell_count)
# - Because the Day 4 subset has one culture level, no
#   culture random effect is included.
# - Main reference level for relative effects is 117 lux.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 06b: MODELS EXP1 LIGHT DAY 4\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check setup objects
# ---------------------------------------------------------

required_objects <- c(
  "project_root",
  "dir_data_derived",
  "dir_tables",
  "dir_logs"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Missing setup object(s):\n- ",
      paste(missing_objects, collapse = "\n- "),
      "\nPlease run 01_setup_packages_and_paths.R first."
    ),
    call. = FALSE
  )
}

cat("Setup objects verified.\n\n")

# ---------------------------------------------------------
# 2. Package checks
# ---------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "stringr",
  "tibble",
  "glmmTMB",
  "broom.mixed"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s):\n- ",
      paste(missing_packages, collapse = "\n- "),
      "\nPlease install these package(s) before running this script."
    ),
    call. = FALSE
  )
}

cat("Required packages verified.\n\n")

# ---------------------------------------------------------
# 3. Source helper scripts
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
# 4. Define input and output paths
# ---------------------------------------------------------

file_input_rds <- file.path(
  dir_data_derived,
  "day4",
  "ms4_day4_analysis_derived.rds"
)

file_input_csv <- file.path(
  dir_data_derived,
  "day4",
  "ms4_day4_analysis_derived.csv"
)

dir_day4_tables <- file.path(project_root, "outputs", "tables", "day4")
dir_day4_logs <- file.path(project_root, "outputs", "logs", "day4")
dir_day4_figures <- file.path(project_root, "outputs", "figures", "day4")
dir_day4_models <- file.path(project_root, "outputs", "models", "day4", "exp1_light_day4")

dir.create(dir_day4_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_models, recursive = TRUE, showWarnings = FALSE)

file_dataset_summary <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_dataset_summary.csv"
)

file_model_comparison <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_model_comparison.csv"
)

file_model_terms <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_preferred_model_terms.csv"
)

file_predictions <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_predictions.csv"
)

file_relative_effects <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_relative_effects.csv"
)

file_preferred_model_rds <- file.path(
  dir_day4_models,
  "exp1_light_day4_preferred_model.rds"
)

file_model_comparison_models <- file.path(
  dir_day4_models,
  "exp1_light_day4_model_comparison.csv"
)

file_predictions_models <- file.path(
  dir_day4_models,
  "exp1_light_day4_predictions.csv"
)

file_relative_effects_models <- file.path(
  dir_day4_models,
  "exp1_light_day4_relative_effects.csv"
)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_log <- file.path(
  dir_day4_logs,
  paste0("06b_exp1_light_day4_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 5. Small local helper functions
# ---------------------------------------------------------

safe_n_distinct <- function(x) {
  dplyr::n_distinct(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

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

inverse_logit <- function(x) {
  stats::plogis(x)
}

fit_model_safely <- function(formula, data, family, model_label) {
  tryCatch(
    {
      model <- glmmTMB::glmmTMB(
        formula = formula,
        data = data,
        family = family
      )
      
      tibble::tibble(
        model_label = model_label,
        model_object = list(model),
        model_status = "success",
        error_message = NA_character_
      )
    },
    error = function(e) {
      tibble::tibble(
        model_label = model_label,
        model_object = list(NULL),
        model_status = "failed",
        error_message = conditionMessage(e)
      )
    }
  )
}

extract_aic_row <- function(model, model_label, family_label, model_status, error_message) {
  
  model_is_valid <- inherits(model, "glmmTMB")
  
  if (is.null(model) || !model_is_valid) {
    return(
      tibble::tibble(
        model = model_label,
        family = family_label,
        df = NA_real_,
        AIC = NA_real_,
        delta_aic = NA_real_,
        model_status = ifelse(model_status == "success", "failed", model_status),
        error_message = dplyr::case_when(
          !is.na(error_message) ~ error_message,
          is.null(model) ~ "Model object is NULL.",
          !model_is_valid ~ paste0(
            "Object passed to AIC extraction is not a glmmTMB model. Class: ",
            paste(class(model), collapse = ", ")
          ),
          TRUE ~ NA_character_
        )
      )
    )
  }
  
  tibble::tibble(
    model = model_label,
    family = family_label,
    df = attr(stats::logLik(model), "df"),
    AIC = stats::AIC(model),
    delta_aic = NA_real_,
    model_status = model_status,
    error_message = error_message
  )
}

predict_response_with_ci <- function(model, newdata, conf_level = 0.95) {
  pred_link <- stats::predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )
  
  z_value <- stats::qnorm(1 - ((1 - conf_level) / 2))
  
  out <- newdata |>
    dplyr::mutate(
      fit_link = pred_link$fit,
      se_link = pred_link$se.fit,
      conf_low_link = fit_link - z_value * se_link,
      conf_high_link = fit_link + z_value * se_link,
      fit_response = inverse_logit(fit_link),
      conf_low_response = inverse_logit(conf_low_link),
      conf_high_response = inverse_logit(conf_high_link)
    )
  
  out
}

save_plot_multi_local <- function(plot, file_stem, width = 7, height = 5.5, dpi = 600) {
  file_pdf <- file.path(dir_day4_figures, paste0(file_stem, ".pdf"))
  file_png <- file.path(dir_day4_figures, paste0(file_stem, ".png"))
  file_tiff <- file.path(dir_day4_figures, paste0(file_stem, ".tiff"))
  
  ggplot2::ggsave(
    filename = file_pdf,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file_png,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file_tiff,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  
  tibble::tibble(
    figure = file_stem,
    file_pdf = file_pdf,
    file_png = file_png,
    file_tiff = file_tiff
  )
}

# ---------------------------------------------------------
# 6. Load Day 4 dataset
# ---------------------------------------------------------

if (file.exists(file_input_rds)) {
  dat <- readRDS(file_input_rds)
  input_source_used <- file_input_rds
} else if (file.exists(file_input_csv)) {
  dat <- readr::read_csv(file_input_csv, show_col_types = FALSE)
  input_source_used <- file_input_csv
} else {
  stop(
    paste0(
      "No Day 4 derived dataset found.\nExpected one of:\n- ",
      file_input_rds,
      "\n- ",
      file_input_csv,
      "\nPlease run 04b_derive_day4_analysis_dataset.R first."
    ),
    call. = FALSE
  )
}

cat("Loaded Day 4 dataset from:\n")
cat(input_source_used, "\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

# ---------------------------------------------------------
# 7. Check required columns
# ---------------------------------------------------------

required_columns <- c(
  "experiment_num",
  "days_from_start",
  "mobile_cell_count",
  "stationary_cell_count",
  "total_cells",
  "motility_ratio",
  "lux_exposure",
  "culture",
  "well"
)

missing_columns <- required_columns[
  !required_columns %in% names(dat)
]

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required column(s):\n- ",
      paste(missing_columns, collapse = "\n- "),
      "\nPlease review the Day 4 derived dataset."
    ),
    call. = FALSE
  )
}

cat("Required columns verified.\n\n")

# ---------------------------------------------------------
# 8. Restrict to Experiment 1: light-only, Day 4
# ---------------------------------------------------------

exp1_dat <- dat |>
  dplyr::mutate(
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    lux_exposure = suppressWarnings(as.numeric(lux_exposure)),
    mobile_cell_count = as.integer(mobile_cell_count),
    stationary_cell_count = as.integer(stationary_cell_count),
    total_cells = mobile_cell_count + stationary_cell_count,
    motility_ratio = dplyr::if_else(
      total_cells > 0,
      mobile_cell_count / total_cells,
      NA_real_
    )
  ) |>
  dplyr::filter(
    experiment_num == 9.2,
    days_from_start == 4,
    !is.na(lux_exposure),
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    total_cells > 0
  )

if (nrow(exp1_dat) == 0) {
  stop("No valid Day 4 Experiment 1 rows found.", call. = FALSE)
}

expected_lux_levels <- c(0, 4, 70, 117)

exp1_dat <- exp1_dat |>
  dplyr::mutate(
    lux_exposure_f = factor(
      lux_exposure,
      levels = c(117, 70, 4, 0),
      labels = c("117 lux", "70 lux", "4 lux", "0 lux")
    ),
    lux_plot_f = factor(
      lux_exposure,
      levels = c(0, 4, 70, 117),
      labels = c("0 lux", "4 lux", "70 lux", "117 lux")
    ),
    culture = as.factor(culture),
    well = as.factor(well),
    experiment_label = "Experiment 1: Light-only gradient",
    reference_condition = "117 lux"
  )

if (any(is.na(exp1_dat$lux_exposure_f))) {
  warning(
    "Some Experiment 1 rows contain lux values outside expected levels: ",
    paste(sort(unique(exp1_dat$lux_exposure)), collapse = ", "),
    call. = FALSE
  )
}

cat("Experiment 1 Day 4 subset created.\n")
cat("Rows:", nrow(exp1_dat), "\n")
cat("Lux levels:", paste(sort(unique(exp1_dat$lux_exposure)), collapse = ", "), "\n")
cat("Cultures:", safe_n_distinct(exp1_dat$culture), "\n")
cat("Wells:", safe_n_distinct(exp1_dat$well), "\n\n")

# ---------------------------------------------------------
# 9. Dataset summary
# ---------------------------------------------------------

has_video_file <- "video_file" %in% names(exp1_dat)

dataset_summary <- exp1_dat |>
  dplyr::group_by(lux_exposure, lux_plot_f) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_cultures = safe_n_distinct(culture),
    n_wells = safe_n_distinct(well),
    n_videos = if (has_video_file) safe_n_distinct(video_file) else NA_integer_,
    total_mobile_cells = sum(mobile_cell_count, na.rm = TRUE),
    total_stationary_cells = sum(stationary_cell_count, na.rm = TRUE),
    total_cells = sum(total_cells, na.rm = TRUE),
    mean_motility_ratio = safe_mean(motility_ratio),
    sd_motility_ratio = safe_sd(motility_ratio),
    min_motility_ratio = safe_min(motility_ratio),
    max_motility_ratio = safe_max(motility_ratio),
    .groups = "drop"
  ) |>
  dplyr::arrange(lux_exposure)

readr::write_csv(dataset_summary, file_dataset_summary)

cat("Dataset summary written to:\n")
cat(file_dataset_summary, "\n\n")

print(dataset_summary)

# ---------------------------------------------------------
# 10. Fit candidate models
# ---------------------------------------------------------

formula_factor <- cbind(mobile_cell_count, stationary_cell_count) ~ lux_exposure_f
formula_null <- cbind(mobile_cell_count, stationary_cell_count) ~ 1

fit_factor_betabinomial <- fit_model_safely(
  formula = formula_factor,
  data = exp1_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "factor_betabinomial"
)

fit_null_betabinomial <- fit_model_safely(
  formula = formula_null,
  data = exp1_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "null_betabinomial"
)

fit_factor_binomial <- fit_model_safely(
  formula = formula_factor,
  data = exp1_dat,
  family = stats::binomial(link = "logit"),
  model_label = "factor_binomial"
)

fit_null_binomial <- fit_model_safely(
  formula = formula_null,
  data = exp1_dat,
  family = stats::binomial(link = "logit"),
  model_label = "null_binomial"
)

fit_rows <- dplyr::bind_rows(
  fit_factor_betabinomial,
  fit_null_betabinomial,
  fit_factor_binomial,
  fit_null_binomial
)

# ---------------------------------------------------------
# 10b. Compare candidate models safely
# ---------------------------------------------------------

candidate_models <- list(
  factor_betabinomial = list(
    object = fit_factor_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_factor_betabinomial$model_status[[1]],
    error = fit_factor_betabinomial$error_message[[1]]
  ),
  null_betabinomial = list(
    object = fit_null_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_null_betabinomial$model_status[[1]],
    error = fit_null_betabinomial$error_message[[1]]
  ),
  factor_binomial = list(
    object = fit_factor_binomial$model_object[[1]],
    family = "binomial",
    status = fit_factor_binomial$model_status[[1]],
    error = fit_factor_binomial$error_message[[1]]
  ),
  null_binomial = list(
    object = fit_null_binomial$model_object[[1]],
    family = "binomial",
    status = fit_null_binomial$model_status[[1]],
    error = fit_null_binomial$error_message[[1]]
  )
)

make_aic_row <- function(model_name, model_info) {
  
  model_object <- model_info$object
  
  if (!inherits(model_object, "glmmTMB")) {
    return(
      tibble::tibble(
        model = model_name,
        family = model_info$family,
        df = NA_real_,
        AIC = NA_real_,
        delta_aic = NA_real_,
        model_status = "failed",
        error_message = paste0(
          "No valid glmmTMB model object available. Object class was: ",
          paste(class(model_object), collapse = ", "),
          ". Original model status: ",
          model_info$status,
          ". Original error: ",
          model_info$error
        )
      )
    )
  }
  
  tibble::tibble(
    model = model_name,
    family = model_info$family,
    df = attr(logLik(model_object), "df"),
    AIC = AIC(model_object),
    delta_aic = NA_real_,
    model_status = "success",
    error_message = NA_character_
  )
}

model_comparison <- dplyr::bind_rows(
  lapply(names(candidate_models), function(nm) {
    make_aic_row(nm, candidate_models[[nm]])
  })
) |>
  dplyr::arrange(AIC) |>
  dplyr::mutate(
    delta_aic = AIC - min(AIC, na.rm = TRUE),
    experiment = "Experiment 1",
    experiment_num = 9.2,
    analysis_branch = "day4"
  )

readr::write_csv(model_comparison, file_model_comparison)
readr::write_csv(model_comparison, file_model_comparison_models)

cat("Model comparison written to:\n")
cat(file_model_comparison, "\n\n")

print(model_comparison)

# ---------------------------------------------------------
# 11. Select preferred model
# ---------------------------------------------------------

preferred_model_name <- model_comparison |>
  dplyr::filter(model_status == "success", !is.na(AIC)) |>
  dplyr::arrange(AIC) |>
  dplyr::slice(1) |>
  dplyr::pull(model)

preferred_model <- switch(
  preferred_model_name,
  factor_betabinomial = fit_factor_betabinomial$model_object[[1]],
  null_betabinomial = fit_null_betabinomial$model_object[[1]],
  factor_binomial = fit_factor_binomial$model_object[[1]],
  null_binomial = fit_null_binomial$model_object[[1]],
  NULL
)

if (is.null(preferred_model)) {
  stop("No preferred model could be selected.", call. = FALSE)
}

saveRDS(preferred_model, file_preferred_model_rds)

cat("Preferred model selected:\n")
cat(preferred_model_name, "\n")
cat("Preferred model saved to:\n")
cat(file_preferred_model_rds, "\n\n")

# ---------------------------------------------------------
# 12. Preferred model fixed-effect table
# ---------------------------------------------------------

preferred_terms <- tryCatch(
  {
    broom.mixed::tidy(
      preferred_model,
      effects = "fixed",
      conf.int = TRUE,
      conf.level = 0.95
    ) |>
      dplyr::mutate(
        model = preferred_model_name,
        experiment = "Experiment 1",
        experiment_num = 9.2,
        analysis_branch = "day4"
      )
  },
  error = function(e) {
    tibble::tibble(
      model = preferred_model_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      error_message = conditionMessage(e),
      experiment = "Experiment 1",
      experiment_num = 9.2,
      analysis_branch = "day4"
    )
  }
)

readr::write_csv(preferred_terms, file_model_terms)

cat("Preferred model terms written to:\n")
cat(file_model_terms, "\n\n")

print(preferred_terms)

# ---------------------------------------------------------
# 13. Prediction grid
# ---------------------------------------------------------

prediction_grid <- tibble::tibble(
  lux_exposure = c(117, 70, 4, 0),
  lux_exposure_f = factor(
    c("117 lux", "70 lux", "4 lux", "0 lux"),
    levels = levels(exp1_dat$lux_exposure_f)
  ),
  lux_plot_f = factor(
    c("117 lux", "70 lux", "4 lux", "0 lux"),
    levels = levels(exp1_dat$lux_plot_f)
  )
)

# If the preferred model is null, prediction still works;
# unused columns in newdata are harmless.
predictions <- predict_response_with_ci(
  model = preferred_model,
  newdata = prediction_grid,
  conf_level = 0.95
) |>
  dplyr::mutate(
    experiment = "Experiment 1",
    experiment_num = 9.2,
    analysis_branch = "day4",
    model = preferred_model_name,
    response = "motile_fraction",
    reference_condition = "117 lux"
  ) |>
  dplyr::arrange(lux_exposure)

readr::write_csv(predictions, file_predictions)
readr::write_csv(predictions, file_predictions_models)

cat("Predictions written to:\n")
cat(file_predictions, "\n\n")

print(predictions)

# ---------------------------------------------------------
# 14. Relative effects versus 117 lux
# ---------------------------------------------------------

reference_prediction <- predictions |>
  dplyr::filter(lux_exposure == 117) |>
  dplyr::slice(1)

if (nrow(reference_prediction) != 1) {
  stop("Could not identify the 117 lux reference prediction.", call. = FALSE)
}

relative_effects <- predictions |>
  dplyr::filter(lux_exposure != 117) |>
  dplyr::mutate(
    contrast_label = paste0(lux_plot_f, " vs 117 lux"),
    treatment_condition = as.character(lux_plot_f),
    reference_condition = "117 lux",
    fit_control = reference_prediction$fit_response,
    fit_treatment = fit_response,
    percent_change = 100 * ((fit_treatment - fit_control) / fit_control),
    conf_low_percent = 100 * ((conf_low_response - fit_control) / fit_control),
    conf_high_percent = 100 * ((conf_high_response - fit_control) / fit_control),
    absolute_change = fit_treatment - fit_control,
    experiment = "Experiment 1",
    experiment_num = 9.2,
    experiment_label = "Light-only gradient",
    analysis_branch = "day4",
    model = preferred_model_name,
    response = "motile_fraction"
  ) |>
  dplyr::select(
    experiment,
    experiment_num,
    experiment_label,
    analysis_branch,
    model,
    response,
    contrast_label,
    reference_condition,
    treatment_condition,
    fit_control,
    fit_treatment,
    percent_change,
    conf_low_percent,
    conf_high_percent,
    absolute_change,
    fit_response,
    conf_low_response,
    conf_high_response
  )

readr::write_csv(relative_effects, file_relative_effects)
readr::write_csv(relative_effects, file_relative_effects_models)

cat("Relative effects written to:\n")
cat(file_relative_effects, "\n\n")

print(relative_effects)

# ---------------------------------------------------------
# 15. Figure 4: Experiment 1 Day 4 predictions
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

plot_raw <- exp1_dat |>
  dplyr::mutate(
    lux_plot_f = factor(
      lux_plot_f,
      levels = c("0 lux", "4 lux", "70 lux", "117 lux")
    )
  )

plot_pred <- predictions |>
  dplyr::mutate(
    lux_plot_f = factor(
      lux_plot_f,
      levels = c("0 lux", "4 lux", "70 lux", "117 lux")
    )
  )

p_exp1 <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data = plot_raw,
    ggplot2::aes(
      x = lux_plot_f,
      y = motility_ratio
    ),
    width = 0.08,
    height = 0,
    alpha = 0.35,
    size = 1.8
  ) +
  ggplot2::geom_errorbar(
    data = plot_pred,
    ggplot2::aes(
      x = lux_plot_f,
      ymin = conf_low_response,
      ymax = conf_high_response
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    data = plot_pred,
    ggplot2::aes(
      x = lux_plot_f,
      y = fit_response
    ),
    size = 2.8
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Experiment 1: Light-only gradient",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Irradiance treatment",
    y = "Motile fraction",
    caption = paste(
      "Points show raw well-level observations.",
      "Black points and error bars show model-predicted motile fraction with 95% confidence intervals.",
      "Reference for proportional effects: 117 lux."
    )
  )

figure_manifest <- save_plot_multi_local(
  plot = p_exp1,
  file_stem = "Fig4_exp1_light_day4_model_predictions",
  width = 7,
  height = 5.5,
  dpi = 600
)

cat("Figure saved:\n")
print(figure_manifest)

# ---------------------------------------------------------
# 16. Write log
# ---------------------------------------------------------

sink(file_log)
cat("Experiment 1 Day 4 model log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input source:\n")
cat(input_source_used, "\n\n")
cat("Rows in Experiment 1 Day 4 dataset:", nrow(exp1_dat), "\n")
cat("Lux levels:\n")
print(sort(unique(exp1_dat$lux_exposure)))
cat("\nCultures:", safe_n_distinct(exp1_dat$culture), "\n")
cat("Wells:", safe_n_distinct(exp1_dat$well), "\n\n")
cat("Model comparison:\n")
print(model_comparison)
cat("\nPreferred model:\n")
cat(preferred_model_name, "\n\n")
cat("Dataset summary path:\n")
cat(file_dataset_summary, "\n\n")
cat("Prediction path:\n")
cat(file_predictions, "\n\n")
cat("Relative effects path:\n")
cat(file_relative_effects, "\n\n")
cat("Figure files:\n")
print(figure_manifest)
sink()

cat("Log written to:\n")
cat(file_log, "\n\n")

cat("SCRIPT 06b COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")