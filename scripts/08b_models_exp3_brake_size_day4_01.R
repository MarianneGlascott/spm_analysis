# =========================================================
# Script title: 08b_models_exp3_brake_size_day4_01.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4 / Article 4
# Purpose: Fit Day 4-only Experiment 3 brake-wear size/
#          fraction models for kelp zoospore motility using
#          primary count data.
#
# Inputs:
# - data_derived/day4/ms4_day4_analysis_derived.rds or .csv
#
# Outputs:
# - outputs/models/day4/exp3_brake_size_day4/*.rds
# - outputs/models/day4/exp3_brake_size_day4/*.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_dataset_summary.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_model_comparison.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_preferred_model_terms.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_predictions.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_relative_effects.csv
# - outputs/figures/day4/Fig6_exp3_brake_size_day4_model_predictions.{pdf,png,tiff}
# - outputs/logs/day4/08b_exp3_brake_size_day4_log_*.txt
#
# Notes:
# - This is the Day 4 analysis branch.
# - This script does not modify the full time-course analysis.
# - Primary response:
#   cbind(mobile_cell_count, stationary_cell_count)
# - Because the Day 4 subset has one culture level, no
#   culture random effect is included.
# - Primary exposure variable:
#   log10_bw_ug_l_plus1 = log10(bw_ug_l + 1)
# - Experiment 3 is interpreted as a brake-wear size/fraction
#   comparison under a brake-wear mass-loading gradient.
# - Main Figure 8 reference condition:
#   BWC at 4.5 ug/L.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 08b: MODELS EXP3 BRAKE SIZE DAY 4\n")
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
  "tidyr",
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
# 3. Source helper scripts if available
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
  } else {
    cat("Helper not found, continuing without:", hf, "\n")
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
dir_day4_models <- file.path(
  project_root,
  "outputs",
  "models",
  "day4",
  "exp3_brake_size_day4"
)

dir.create(dir_day4_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_models, recursive = TRUE, showWarnings = FALSE)

file_dataset_summary <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_dataset_summary.csv"
)

file_model_comparison <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_model_comparison.csv"
)

file_model_terms <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_preferred_model_terms.csv"
)

file_predictions <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_predictions.csv"
)

file_relative_effects <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_relative_effects.csv"
)

file_preferred_model_rds <- file.path(
  dir_day4_models,
  "exp3_brake_size_day4_preferred_model.rds"
)

file_model_comparison_models <- file.path(
  dir_day4_models,
  "exp3_brake_size_day4_model_comparison.csv"
)

file_predictions_models <- file.path(
  dir_day4_models,
  "exp3_brake_size_day4_predictions.csv"
)

file_relative_effects_models <- file.path(
  dir_day4_models,
  "exp3_brake_size_day4_relative_effects.csv"
)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_log <- file.path(
  dir_day4_logs,
  paste0("08b_exp3_brake_size_day4_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 5. Small local helper functions
# ---------------------------------------------------------

standardise_upper_trim <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
}

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

save_plot_multi_local <- function(plot, file_stem, width = 7.8, height = 5.8, dpi = 600) {
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

# bw_ug_l is required for this script. If not present, stop explicitly.
if (!"bw_ug_l" %in% names(dat)) {
  stop(
    "Column bw_ug_l is required for Experiment 3 brake-wear size modelling but is missing.",
    call. = FALSE
  )
}

cat("Required columns verified.\n\n")

# ---------------------------------------------------------
# 8. Restrict to Experiment 3: brake-wear size, Day 4
# ---------------------------------------------------------

exp3_dat <- dat |>
  dplyr::mutate(
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    bw_ug_l = suppressWarnings(as.numeric(bw_ug_l)),
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
    experiment_num == 11.2,
    days_from_start == 4,
    !is.na(bw_ug_l),
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    total_cells > 0
  )

if (nrow(exp3_dat) == 0) {
  stop("No valid Day 4 Experiment 3 rows found.", call. = FALSE)
}

# Resolve size/fraction identity from available upstream columns.
size_source_col <- dplyr::case_when(
  "size_class_plot" %in% names(exp3_dat) ~ "size_class_plot",
  "size_class" %in% names(exp3_dat) ~ "size_class",
  "brake_size_treatment" %in% names(exp3_dat) ~ "brake_size_treatment",
  "toxin_exposure" %in% names(exp3_dat) ~ "toxin_exposure",
  TRUE ~ NA_character_
)

if (is.na(size_source_col)) {
  stop(
    "No brake-wear size/fraction column found. Expected one of size_class_plot, size_class, brake_size_treatment, or toxin_exposure.",
    call. = FALSE
  )
}

exp3_dat$size_source <- exp3_dat[[size_source_col]]

exp3_dat <- exp3_dat |>
  dplyr::mutate(
    particle_size_class = standardise_upper_trim(size_source),
    particle_size_class = dplyr::case_when(
      particle_size_class %in% c("BWC", "COARSE", "BRAKE WEAR COARSE", "BRAKE_WEAR_COARSE", "BRAKE COARSE") ~ "BWC",
      particle_size_class %in% c("BWF", "FINE", "BRAKE WEAR FINE", "BRAKE_WEAR_FINE", "BRAKE FINE") ~ "BWF",
      TRUE ~ particle_size_class
    ),
    log10_bw_ug_l_plus1 = log10(bw_ug_l + 1),
    bw_load_f = factor(
      bw_ug_l,
      levels = c(4.5, 45, 450),
      labels = c("4.5 ug/L", "45 ug/L", "450 ug/L")
    ),
    culture = as.factor(culture),
    well = as.factor(well),
    experiment_label = "Experiment 3: Brake-wear size comparison",
    analysis_branch = "day4"
  ) |>
  dplyr::filter(
    particle_size_class %in% c("BWC", "BWF")
  ) |>
  dplyr::mutate(
    particle_size_class = factor(
      particle_size_class,
      levels = c("BWC", "BWF")
    ),
    particle_size_plot = factor(
      particle_size_class,
      levels = c("BWC", "BWF"),
      labels = c("Coarse", "Fine")
    )
  )

if (nrow(exp3_dat) == 0) {
  stop("No valid BWC/BWF rows remain for Experiment 3.", call. = FALSE)
}

if (safe_n_distinct(exp3_dat$particle_size_class) < 2) {
  stop(
    "Experiment 3 requires both BWC and BWF levels for the size/fraction comparison.",
    call. = FALSE
  )
}

cat("Experiment 3 Day 4 subset created.\n")
cat("Rows:", nrow(exp3_dat), "\n")
cat("Size/fraction levels:", paste(levels(droplevels(exp3_dat$particle_size_class)), collapse = ", "), "\n")
cat("Brake-wear loads (ug/L):", paste(sort(unique(exp3_dat$bw_ug_l)), collapse = ", "), "\n")
cat("Culture levels:", safe_n_distinct(exp3_dat$culture), "\n")
cat("Wells:", safe_n_distinct(exp3_dat$well), "\n")
cat("Size/fraction source column:", size_source_col, "\n\n")

# ---------------------------------------------------------
# 9. Dataset summary
# ---------------------------------------------------------

has_video_file <- "video_file" %in% names(exp3_dat)

dataset_summary <- exp3_dat |>
  dplyr::group_by(particle_size_class, particle_size_plot, bw_ug_l, bw_load_f) |>
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
  dplyr::arrange(particle_size_class, bw_ug_l)

readr::write_csv(dataset_summary, file_dataset_summary)

cat("Dataset summary written to:\n")
cat(file_dataset_summary, "\n\n")

print(dataset_summary)

# ---------------------------------------------------------
# 10. Fit candidate models
# ---------------------------------------------------------

formula_interaction <- cbind(mobile_cell_count, stationary_cell_count) ~
  particle_size_class * log10_bw_ug_l_plus1

formula_additive <- cbind(mobile_cell_count, stationary_cell_count) ~
  particle_size_class + log10_bw_ug_l_plus1

formula_size_only <- cbind(mobile_cell_count, stationary_cell_count) ~
  particle_size_class

formula_concentration_only <- cbind(mobile_cell_count, stationary_cell_count) ~
  log10_bw_ug_l_plus1

formula_null <- cbind(mobile_cell_count, stationary_cell_count) ~ 1

fit_interaction_betabinomial <- fit_model_safely(
  formula = formula_interaction,
  data = exp3_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "interaction_betabinomial"
)

fit_additive_betabinomial <- fit_model_safely(
  formula = formula_additive,
  data = exp3_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "additive_betabinomial"
)

fit_size_only_betabinomial <- fit_model_safely(
  formula = formula_size_only,
  data = exp3_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "size_only_betabinomial"
)

fit_concentration_only_betabinomial <- fit_model_safely(
  formula = formula_concentration_only,
  data = exp3_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "concentration_only_betabinomial"
)

fit_null_betabinomial <- fit_model_safely(
  formula = formula_null,
  data = exp3_dat,
  family = glmmTMB::betabinomial(link = "logit"),
  model_label = "null_betabinomial"
)

fit_interaction_binomial <- fit_model_safely(
  formula = formula_interaction,
  data = exp3_dat,
  family = stats::binomial(link = "logit"),
  model_label = "interaction_binomial"
)

fit_additive_binomial <- fit_model_safely(
  formula = formula_additive,
  data = exp3_dat,
  family = stats::binomial(link = "logit"),
  model_label = "additive_binomial"
)

fit_size_only_binomial <- fit_model_safely(
  formula = formula_size_only,
  data = exp3_dat,
  family = stats::binomial(link = "logit"),
  model_label = "size_only_binomial"
)

fit_concentration_only_binomial <- fit_model_safely(
  formula = formula_concentration_only,
  data = exp3_dat,
  family = stats::binomial(link = "logit"),
  model_label = "concentration_only_binomial"
)

fit_null_binomial <- fit_model_safely(
  formula = formula_null,
  data = exp3_dat,
  family = stats::binomial(link = "logit"),
  model_label = "null_binomial"
)

candidate_models <- list(
  interaction_betabinomial = list(
    object = fit_interaction_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_interaction_betabinomial$model_status[[1]],
    error = fit_interaction_betabinomial$error_message[[1]]
  ),
  additive_betabinomial = list(
    object = fit_additive_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_additive_betabinomial$model_status[[1]],
    error = fit_additive_betabinomial$error_message[[1]]
  ),
  size_only_betabinomial = list(
    object = fit_size_only_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_size_only_betabinomial$model_status[[1]],
    error = fit_size_only_betabinomial$error_message[[1]]
  ),
  concentration_only_betabinomial = list(
    object = fit_concentration_only_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_concentration_only_betabinomial$model_status[[1]],
    error = fit_concentration_only_betabinomial$error_message[[1]]
  ),
  null_betabinomial = list(
    object = fit_null_betabinomial$model_object[[1]],
    family = "betabinomial",
    status = fit_null_betabinomial$model_status[[1]],
    error = fit_null_betabinomial$error_message[[1]]
  ),
  interaction_binomial = list(
    object = fit_interaction_binomial$model_object[[1]],
    family = "binomial",
    status = fit_interaction_binomial$model_status[[1]],
    error = fit_interaction_binomial$error_message[[1]]
  ),
  additive_binomial = list(
    object = fit_additive_binomial$model_object[[1]],
    family = "binomial",
    status = fit_additive_binomial$model_status[[1]],
    error = fit_additive_binomial$error_message[[1]]
  ),
  size_only_binomial = list(
    object = fit_size_only_binomial$model_object[[1]],
    family = "binomial",
    status = fit_size_only_binomial$model_status[[1]],
    error = fit_size_only_binomial$error_message[[1]]
  ),
  concentration_only_binomial = list(
    object = fit_concentration_only_binomial$model_object[[1]],
    family = "binomial",
    status = fit_concentration_only_binomial$model_status[[1]],
    error = fit_concentration_only_binomial$error_message[[1]]
  ),
  null_binomial = list(
    object = fit_null_binomial$model_object[[1]],
    family = "binomial",
    status = fit_null_binomial$model_status[[1]],
    error = fit_null_binomial$error_message[[1]]
  )
)

model_comparison <- dplyr::bind_rows(
  lapply(names(candidate_models), function(nm) {
    make_aic_row(nm, candidate_models[[nm]])
  })
) |>
  dplyr::arrange(AIC) |>
  dplyr::mutate(
    delta_aic = AIC - min(AIC, na.rm = TRUE),
    experiment = "Experiment 3",
    experiment_num = 11.2,
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

preferred_model <- candidate_models[[preferred_model_name]]$object

if (is.null(preferred_model) || !inherits(preferred_model, "glmmTMB")) {
  stop("No valid preferred model could be selected.", call. = FALSE)
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
        experiment = "Experiment 3",
        experiment_num = 11.2,
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
      experiment = "Experiment 3",
      experiment_num = 11.2,
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

observed_bw_loads <- sort(unique(exp3_dat$bw_ug_l))

expected_bw_loads <- c(4.5, 45, 450)

if (all(expected_bw_loads %in% observed_bw_loads)) {
  prediction_bw_loads <- expected_bw_loads
} else {
  prediction_bw_loads <- observed_bw_loads
  warning(
    paste0(
      "Expected brake-wear loads 4.5, 45, 450 ug/L were not all observed. ",
      "Using observed loads instead: ",
      paste(prediction_bw_loads, collapse = ", ")
    ),
    call. = FALSE
  )
}

prediction_grid <- tidyr::expand_grid(
  particle_size_class = factor(
    c("BWC", "BWF"),
    levels = levels(exp3_dat$particle_size_class)
  ),
  bw_ug_l = prediction_bw_loads
) |>
  dplyr::mutate(
    log10_bw_ug_l_plus1 = log10(bw_ug_l + 1),
    bw_load_f = factor(
      bw_ug_l,
      levels = prediction_bw_loads,
      labels = paste0(prediction_bw_loads, " ug/L")
    ),
    particle_size_plot = factor(
      particle_size_class,
      levels = c("BWC", "BWF"),
      labels = c("Coarse", "Fine")
    )
  )

predictions <- predict_response_with_ci(
  model = preferred_model,
  newdata = prediction_grid,
  conf_level = 0.95
) |>
  dplyr::mutate(
    experiment = "Experiment 3",
    experiment_num = 11.2,
    experiment_label = "Brake-wear size comparison",
    analysis_branch = "day4",
    model = preferred_model_name,
    response = "motile_fraction",
    reference_condition = "BWC at 4.5 ug/L",
    treatment_condition = paste0(
      as.character(particle_size_plot),
      " at ",
      bw_ug_l,
      " ug/L"
    )
  ) |>
  dplyr::arrange(particle_size_class, bw_ug_l)

readr::write_csv(predictions, file_predictions)
readr::write_csv(predictions, file_predictions_models)

cat("Predictions written to:\n")
cat(file_predictions, "\n\n")

print(predictions)

# ---------------------------------------------------------
# 14. Relative effects versus BWC at 4.5 ug/L
# ---------------------------------------------------------

reference_bw_load <- 4.5

if (!reference_bw_load %in% predictions$bw_ug_l) {
  reference_bw_load <- min(predictions$bw_ug_l, na.rm = TRUE)
  warning(
    paste0(
      "BWC at 4.5 ug/L was not available in predictions. ",
      "Using BWC at lowest predicted load instead: ",
      reference_bw_load,
      " ug/L."
    ),
    call. = FALSE
  )
}

reference_prediction <- predictions |>
  dplyr::filter(
    particle_size_class == "BWC",
    bw_ug_l == reference_bw_load
  ) |>
  dplyr::slice(1)

if (nrow(reference_prediction) != 1) {
  stop("Could not identify the BWC low-load reference prediction.", call. = FALSE)
}

reference_condition_label <- paste0("BWC at ", reference_bw_load, " ug/L")

relative_effects <- predictions |>
  dplyr::filter(!(particle_size_class == "BWC" & bw_ug_l == reference_bw_load)) |>
  dplyr::mutate(
    contrast_label = paste0(treatment_condition, " vs ", reference_condition_label),
    reference_condition = reference_condition_label,
    fit_control = reference_prediction$fit_response,
    fit_treatment = fit_response,
    percent_change = 100 * ((fit_treatment - fit_control) / fit_control),
    conf_low_percent = 100 * ((conf_low_response - fit_control) / fit_control),
    conf_high_percent = 100 * ((conf_high_response - fit_control) / fit_control),
    absolute_change = fit_treatment - fit_control,
    experiment = "Experiment 3",
    experiment_num = 11.2,
    experiment_label = "Brake-wear size comparison",
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
    particle_size_class,
    particle_size_plot,
    bw_ug_l,
    bw_load_f,
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
# 15. Figure 6: Experiment 3 Day 4 predictions
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

size_palette_local <- c(
  "Coarse" = "#6E6E6E",
  "Fine" = "#2F2F2F"
)

plot_raw <- exp3_dat |>
  dplyr::mutate(
    bw_load_f = factor(
      bw_ug_l,
      levels = prediction_bw_loads,
      labels = paste0(prediction_bw_loads, " ug/L")
    ),
    particle_size_plot = factor(
      particle_size_plot,
      levels = c("Coarse", "Fine")
    )
  )

plot_pred <- predictions |>
  dplyr::mutate(
    bw_load_f = factor(
      bw_ug_l,
      levels = prediction_bw_loads,
      labels = paste0(prediction_bw_loads, " ug/L")
    ),
    particle_size_plot = factor(
      particle_size_plot,
      levels = c("Coarse", "Fine")
    )
  )

p_exp3 <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data = plot_raw,
    ggplot2::aes(
      x = bw_load_f,
      y = motility_ratio,
      colour = particle_size_plot
    ),
    width = 0.10,
    height = 0,
    alpha = 0.35,
    size = 1.8
  ) +
  ggplot2::geom_errorbar(
    data = plot_pred,
    ggplot2::aes(
      x = bw_load_f,
      ymin = conf_low_response,
      ymax = conf_high_response,
      colour = particle_size_plot
    ),
    width = 0.15,
    linewidth = 0.6,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::geom_point(
    data = plot_pred,
    ggplot2::aes(
      x = bw_load_f,
      y = fit_response,
      colour = particle_size_plot
    ),
    size = 2.8,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::scale_colour_manual(values = size_palette_local, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Experiment 3: Brake-wear size comparison",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Brake-wear load",
    y = "Motile fraction",
    colour = "Brake-wear fraction",
    caption = paste(
      "Points show raw well-level observations.",
      "Model predictions are shown for coarse and fine brake-wear fractions across mass loads.",
      "Reference for proportional effects: BWC at 4.5 ug/L."
    )
  )

figure_manifest <- save_plot_multi_local(
  plot = p_exp3,
  file_stem = "Fig6_exp3_brake_size_day4_model_predictions",
  width = 7.8,
  height = 5.8,
  dpi = 600
)

cat("Figure saved:\n")
print(figure_manifest)

# ---------------------------------------------------------
# 16. Write log
# ---------------------------------------------------------

sink(file_log)
cat("Experiment 3 Day 4 model log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input source:\n")
cat(input_source_used, "\n\n")
cat("Rows in Experiment 3 Day 4 dataset:", nrow(exp3_dat), "\n")
cat("Size/fraction source column:", size_source_col, "\n\n")
cat("Size/fraction levels:\n")
print(table(exp3_dat$particle_size_class, useNA = "ifany"))
cat("\nBrake-wear loads:\n")
print(sort(unique(exp3_dat$bw_ug_l)))
cat("\nCultures:", safe_n_distinct(exp3_dat$culture), "\n")
cat("Wells:", safe_n_distinct(exp3_dat$well), "\n\n")
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

cat("SCRIPT 08b COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")