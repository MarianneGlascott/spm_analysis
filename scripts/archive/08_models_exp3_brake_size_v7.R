# =========================================================
# Script title: 08_models_exp3_brake_size_v7.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Fit and evaluate Experiment 3 brake-wear size
#          models for motile fraction using primary cell
#          counts, explicitly incorporating brake-wear
#          concentration (bw_ug_l) as the primary exposure
#          gradient, while adapting automatically to the
#          observed design (e.g. single sampled day or
#          single culture level).
# Inputs:
# - data_derived/ms4_analysis_derived.csv or .rds
# Outputs:
# - outputs/models/exp3_brake_size/*.rds
# - outputs/models/exp3_brake_size/*_fixed_effects.csv
# - outputs/models/exp3_brake_size/*_metadata.csv
# - outputs/models/exp3_brake_size/*_diagnostics_summary.csv
# - outputs/models/exp3_brake_size/*_model_comparison.csv
# - outputs/models/exp3_brake_size/*_predictions.csv
# - outputs/tables/08_exp3_brake_size_dataset_summary.csv
# - outputs/tables/08_exp3_brake_size_model_comparison.csv
# - outputs/tables/08_exp3_brake_size_model_summary_table.csv
# - outputs/figures/models_exp3/Fig4_brake_size_model_predictions_*.{pdf,png,tiff}
# - outputs/logs/08_models_exp3_brake_size_log_*.txt
# Date created: 28 February 2026
# Last updated: 30 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 04_derive_variables.R before this script.
# - Primary fitted response:
#   cbind(mobile_cell_count, stationary_cell_count)
# - Manuscript 4 block must be restricted to:
#   experiment_num %in% c(8.2, 9.2, 10.2, 11.2)
# - Focal modelling subset for this script:
#   experiment_num == 11.2
# - Experiment 3 design:
#   Brake wear coarse (BWC) vs brake wear fine (BWF),
#   under controlled mass loading using bw_ug_l.
# - bw_ug_l is the primary exposure variable in this
#   experiment and must be modelled on the log10(bw_ug_l + 1)
#   scale.
# - NTU is not the primary predictor in this experiment.
# - If only one observed day is present, days_from_start is
#   retained as metadata/plot label only and is not fitted.
# - If only one culture level is present, no random effect
#   is included for culture.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 08: MODELS EXP3 BRAKE SIZE\n")
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

standardise_upper_trim <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
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
  
  caption_file <- file.path(
    fig_dir,
    paste0(figure_name, "_caption.md")
  )
  
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

# ---------------------------------------------------------
# 4. Define input and output paths
# ---------------------------------------------------------

file_input_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")

exp3_model_subdir <- "exp3_brake_size"
exp3_figure_subdir <- "models_exp3"

dir_model_output <- file.path(project_root, "outputs", "models", exp3_model_subdir)

file_dataset_summary <- file.path(dir_tables, "08_exp3_brake_size_dataset_summary.csv")
file_model_comparison_table <- file.path(dir_tables, "08_exp3_brake_size_model_comparison.csv")
file_model_summary_table <- file.path(dir_tables, "08_exp3_brake_size_model_summary_table.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_model_log <- file.path(
  dir_logs,
  paste0("08_models_exp3_brake_size_log_", timestamp_now, ".txt")
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
  "culture",
  "bw_ug_l",
  "log10_bw_ug_l_plus1"
)

missing_required_columns <- required_columns[!required_columns %in% names(dat)]

if (length(missing_required_columns) > 0) {
  stop(
    paste0(
      "The following required Experiment 3 column(s) are missing:\n- ",
      paste(missing_required_columns, collapse = "\n- "),
      "\nPlease review upstream scripts."
    ),
    call. = FALSE
  )
}

cat("Required Experiment 3 columns verified.\n\n")

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

exp3_dat <- dat |>
  dplyr::filter(experiment_num == 11.2)

if ("use_exp3_brake_size" %in% names(exp3_dat)) {
  exp3_dat <- exp3_dat |>
    dplyr::filter(isTRUE(use_exp3_brake_size))
}

# Prefer upstream derived plotting variables when present
size_source_col <- dplyr::case_when(
  "size_class_plot" %in% names(exp3_dat) ~ "size_class_plot",
  "size_class" %in% names(exp3_dat) ~ "size_class",
  TRUE ~ NA_character_
)

if (is.na(size_source_col)) {
  stop(
    "No Experiment 3 size-class column was found (expected size_class_plot or size_class).",
    call. = FALSE
  )
}

if (!"days_from_start" %in% names(exp3_dat)) {
  exp3_dat$days_from_start <- NA_integer_
}

exp3_dat <- exp3_dat |>
  dplyr::mutate(
    size_class_source = .data[[size_source_col]],
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    culture = as.factor(culture),
    bw_ug_l = suppressWarnings(as.numeric(bw_ug_l)),
    log10_bw_ug_l_plus1 = suppressWarnings(as.numeric(log10_bw_ug_l_plus1))
  ) |>
  dplyr::filter(
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    !is.na(culture),
    !is.na(size_class_source),
    !is.na(bw_ug_l),
    !is.na(log10_bw_ug_l_plus1)
  )

if (nrow(exp3_dat) == 0) {
  stop("No Experiment 3 rows available after filtering.", call. = FALSE)
}

cat("Experiment 3 focal subset created.\n")
cat("Rows:", nrow(exp3_dat), "\n")
cat("Cultures:", dplyr::n_distinct(exp3_dat$culture, na.rm = TRUE), "\n")
cat(
  "Days:",
  if ("days_from_start" %in% names(exp3_dat)) {
    paste(sort(unique(stats::na.omit(exp3_dat$days_from_start))), collapse = ", ")
  } else {
    "not available"
  },
  "\n\n"
)

# ---------------------------------------------------------
# 8. Standardise size class and derive support vars
# ---------------------------------------------------------

exp3_dat <- exp3_dat |>
  dplyr::mutate(
    particle_size_class = standardise_upper_trim(size_class_source),
    particle_size_class = dplyr::case_when(
      particle_size_class %in% c("COARSE", "BWC", "BRAKE WEAR COARSE", "BRAKE_WEAR_COARSE") ~ "BWC",
      particle_size_class %in% c("FINE", "BWF", "BRAKE WEAR FINE", "BRAKE_WEAR_FINE") ~ "BWF",
      TRUE ~ NA_character_
    ),
    total_cells = mobile_cell_count + stationary_cell_count,
    motility_ratio = dplyr::if_else(
      total_cells > 0,
      mobile_cell_count / total_cells,
      NA_real_
    ),
    day_label = dplyr::if_else(
      !is.na(days_from_start),
      paste0("Day ", days_from_start),
      "Single sampling day"
    ),
    experiment_plot_label = "Experiment 3: Brake size"
  ) |>
  dplyr::filter(
    !is.na(particle_size_class),
    particle_size_class %in% c("BWC", "BWF"),
    total_cells > 0
  ) |>
  dplyr::mutate(
    particle_size_class = factor(
      particle_size_class,
      levels = c("BWC", "BWF")
    )
  )

if (nrow(exp3_dat) == 0) {
  stop("No valid BWC/BWF rows remain after standardising size class.", call. = FALSE)
}

observed_size_levels <- levels(droplevels(exp3_dat$particle_size_class))
observed_days <- if ("days_from_start" %in% names(exp3_dat)) sort(unique(stats::na.omit(exp3_dat$days_from_start))) else integer(0)
n_day_levels <- length(observed_days)
n_culture_levels <- dplyr::n_distinct(exp3_dat$culture, na.rm = TRUE)

has_day_variation <- n_day_levels > 1
has_culture_variation <- n_culture_levels > 1

cat("Experiment 3 modelling support variables derived.\n")
cat("Size classes observed:", paste(observed_size_levels, collapse = ", "), "\n")
cat(
  "Brake-wear concentration range (ug/L):",
  paste0(
    round(min(exp3_dat$bw_ug_l, na.rm = TRUE), 3),
    " to ",
    round(max(exp3_dat$bw_ug_l, na.rm = TRUE), 3)
  ),
  "\n"
)
cat("Day variation available:", has_day_variation, "\n")
cat("Culture variation available:", has_culture_variation, "\n\n")

# ---------------------------------------------------------
# 9. Build dataset summary table
# ---------------------------------------------------------

has_well <- "well" %in% names(exp3_dat)
has_video_file <- "video_file" %in% names(exp3_dat)

group_vars <- c("particle_size_class")
if ("days_from_start" %in% names(exp3_dat) && !all(is.na(exp3_dat$days_from_start))) {
  group_vars <- c("days_from_start", group_vars)
}

exp3_dataset_summary <- exp3_dat |>
  dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
    n_wells = if (has_well) dplyr::n_distinct(well, na.rm = TRUE) else NA_integer_,
    n_videos = if (has_video_file) dplyr::n_distinct(video_file, na.rm = TRUE) else NA_integer_,
    min_bw_ug_l = min(bw_ug_l, na.rm = TRUE),
    max_bw_ug_l = max(bw_ug_l, na.rm = TRUE),
    mean_bw_ug_l = safe_mean(bw_ug_l),
    mean_total_cells = safe_mean(total_cells),
    mean_motility_ratio = safe_mean(motility_ratio),
    sd_motility_ratio = safe_sd(motility_ratio),
    .groups = "drop"
  )

readr::write_csv(exp3_dataset_summary, file_dataset_summary)

cat("Experiment 3 dataset summary written to:\n")
cat(file_dataset_summary, "\n\n")

# ---------------------------------------------------------
# 10. Build model formulas dynamically from observed design
# ---------------------------------------------------------

rhs_with_random <- function(fixed_rhs, include_random_culture) {
  if (isTRUE(include_random_culture)) {
    paste0(fixed_rhs, " + (1 | culture)")
  } else {
    fixed_rhs
  }
}

as_binomial_formula <- function(rhs) {
  stats::as.formula(
    paste0(
      "cbind(mobile_cell_count, stationary_cell_count) ~ ",
      rhs
    )
  )
}

candidate_formula_strings <- list()

if (has_day_variation) {
  candidate_formula_strings$full_primary <- rhs_with_random(
    "particle_size_class * log10_bw_ug_l_plus1 * days_from_start",
    has_culture_variation
  )
  candidate_formula_strings$partial_primary <- rhs_with_random(
    "particle_size_class * log10_bw_ug_l_plus1 + days_from_start",
    has_culture_variation
  )
  candidate_formula_strings$additive_primary <- rhs_with_random(
    "particle_size_class + log10_bw_ug_l_plus1 + days_from_start",
    has_culture_variation
  )
} else {
  candidate_formula_strings$full_primary <- rhs_with_random(
    "particle_size_class * log10_bw_ug_l_plus1",
    has_culture_variation
  )
  candidate_formula_strings$additive_primary <- rhs_with_random(
    "particle_size_class + log10_bw_ug_l_plus1",
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

for (nm in names(candidate_formulas)) {
  fit_results[[nm]] <- fit_check_save_model(
    formula = candidate_formulas[[nm]],
    data = exp3_dat,
    model_name = paste0("exp3_brake_size_", nm),
    model_subdir = exp3_model_subdir,
    save_terms_csv = TRUE,
    save_meta_csv = TRUE,
    save_diagnostics = TRUE,
    run_dharma = TRUE,
    run_performance = TRUE,
    conf.level = 0.95,
    exponentiate = FALSE,
    quiet = FALSE
  )
}

cat("Candidate models fitted.\n\n")

required_models <- names(candidate_formulas)
missing_models <- required_models[!vapply(required_models, function(x) !is.null(fit_results[[x]]), logical(1))]

if (length(missing_models) > 0) {
  stop(
    paste(
      "The following models failed to fit:",
      paste(missing_models, collapse = ", ")
    ),
    call. = FALSE
  )
}

cat("All candidate models successfully fitted.\n\n")

# ---------------------------------------------------------
# 12. Compare candidate models
# ---------------------------------------------------------

if (!exists("compare_models_aic", mode = "function", inherits = TRUE)) {
  stop("compare_models_aic() is not available.", call. = FALSE)
}

comparison_input <- lapply(fit_results, function(x) x$model)

model_comparison <- do.call(compare_models_aic, comparison_input) |>
  dplyr::mutate(
    experiment = "Experiment 3",
    experiment_num = 11.2,
    complexity_rank = dplyr::case_when(
      model == "full_primary" ~ 4L,
      model == "partial_primary" ~ 3L,
      model == "additive_primary" ~ 2L,
      model == "null" ~ 1L,
      TRUE ~ 99L
    )
  ) |>
  dplyr::arrange(delta_aic, complexity_rank)

readr::write_csv(model_comparison, file_model_comparison_table)

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = model_comparison,
    file_stem = "exp3_brake_size_model_comparison",
    subdir = exp3_model_subdir,
    quiet = FALSE
  )
}

cat("Model comparison written to:\n")
cat(file_model_comparison_table, "\n\n")

# ---------------------------------------------------------
# 13. Likelihood ratio comparisons
# ---------------------------------------------------------

lrt_additive_vs_full <- NULL
lrt_null_vs_additive <- NULL
lrt_partial_vs_full <- NULL

if (all(c("additive_primary", "full_primary") %in% names(fit_results))) {
  lrt_additive_vs_full <- tryCatch(
    stats::anova(fit_results$additive_primary$model, fit_results$full_primary$model),
    error = function(e) e
  )
}

if (all(c("null", "additive_primary") %in% names(fit_results))) {
  lrt_null_vs_additive <- tryCatch(
    stats::anova(fit_results$null$model, fit_results$additive_primary$model),
    error = function(e) e
  )
}

if (all(c("partial_primary", "full_primary") %in% names(fit_results))) {
  lrt_partial_vs_full <- tryCatch(
    stats::anova(fit_results$partial_primary$model, fit_results$full_primary$model),
    error = function(e) e
  )
}

# ---------------------------------------------------------
# 14. Select preferred model
# ---------------------------------------------------------

best_delta <- min(model_comparison$delta_aic, na.rm = TRUE)

preferred_model_name <- model_comparison |>
  dplyr::filter(delta_aic == best_delta) |>
  dplyr::arrange(complexity_rank) |>
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

pred_at <- list(
  particle_size_class = levels(droplevels(exp3_dat$particle_size_class))
)

if (has_day_variation) {
  pred_at$days_from_start <- sort(unique(exp3_dat$days_from_start))
}

pred_grid <- make_prediction_grid(
  data = exp3_dat,
  focal_terms = c("log10_bw_ug_l_plus1", "particle_size_class"),
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
    bw_ug_l_backtransformed = (10^log10_bw_ug_l_plus1) - 1,
    day_label = dplyr::case_when(
      "days_from_start" %in% names(.) & !is.na(days_from_start) ~ paste0("Day ", days_from_start),
      TRUE ~ "Day 4"
    ),
    prediction_scale = "log10_bw_ug_l_plus1"
  )

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = pred_out,
    file_stem = "exp3_brake_size_preferred_model_predictions",
    subdir = exp3_model_subdir,
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

preferred_fixed_effects <- build_model_term_summary(
  model = preferred_model,
  model_name = preferred_model_name,
  conf.level = 0.95,
  exponentiate = FALSE
) |>
  dplyr::mutate(
    experiment = "Experiment 3",
    experiment_num = 11.2
  )

readr::write_csv(preferred_fixed_effects, file_model_summary_table)

cat("Preferred model fixed-effect summary written to:\n")
cat(file_model_summary_table, "\n\n")

# ---------------------------------------------------------
# 17. Build Experiment 3 figure
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

if (exists("particle_palette", mode = "function", inherits = TRUE)) {
  exp3_particle_palette <- particle_palette()
} else if (exists("particle_palette", inherits = TRUE)) {
  exp3_particle_palette <- get("particle_palette", inherits = TRUE)
} else {
  exp3_particle_palette <- c(
    "BWC" = "#D55E00",
    "BWF" = "#CC79A7"
  )
}

# Restrict palette to observed levels if broader project palette is present
exp3_particle_palette <- exp3_particle_palette[names(exp3_particle_palette) %in% c("BWC", "BWF")]
if (length(exp3_particle_palette) == 0) {
  exp3_particle_palette <- c(
    "BWC" = "#D55E00",
    "BWF" = "#CC79A7"
  )
}

plot_raw <- exp3_dat |>
  dplyr::mutate(
    day_label = dplyr::case_when(
      !is.na(days_from_start) ~ paste0("Day ", days_from_start),
      TRUE ~ "Day 4"
    )
  )

p_exp3 <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = plot_raw,
    ggplot2::aes(
      x = log10_bw_ug_l_plus1,
      y = motility_ratio,
      colour = particle_size_class
    ),
    alpha = 0.30,
    size = 1.8
  ) +
  ggplot2::geom_ribbon(
    data = pred_out,
    ggplot2::aes(
      x = log10_bw_ug_l_plus1,
      ymin = conf_low_response,
      ymax = conf_high_response,
      fill = particle_size_class
    ),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = pred_out,
    ggplot2::aes(
      x = log10_bw_ug_l_plus1,
      y = fit_response,
      colour = particle_size_class
    ),
    linewidth = 0.8
  ) +
  ggplot2::scale_colour_manual(values = exp3_particle_palette, drop = FALSE) +
  ggplot2::scale_fill_manual(values = exp3_particle_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Brake-wear size comparison under controlled mass loading",
    subtitle = if (has_day_variation) {
      "Model-predicted motile fraction versus log10(bw_ug_l + 1), by size class and day"
    } else {
      "Model-predicted motile fraction versus log10(bw_ug_l + 1), by size class"
    },
    x = "log10(brake-wear concentration [ug/L] + 1)",
    y = "Motile fraction",
    colour = "Particle class",
    fill = "Particle class",
    caption = if (has_day_variation) {
      "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals. Panels are shown by day from start."
    } else {
      "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals. Experiment 3 was sampled on a single day only, so no day facet is shown."
    }
  )

if (has_day_variation) {
  p_exp3 <- p_exp3 + ggplot2::facet_wrap(~ day_label)
}

# ---------------------------------------------------------
# 18. Save Experiment 3 figure
# ---------------------------------------------------------

saved_figures <- character(0)

if (exists("save_figure_both_widths", mode = "function", inherits = TRUE)) {
  saved_figures <- save_figure_both_widths(
    plot = p_exp3,
    figure_name = "Fig4_brake_size_model_predictions",
    subdir = exp3_figure_subdir,
    height = "standard",
    quiet = TRUE
  )
}

caption_text <- if (has_day_variation) {
  "Model-predicted motile fraction as a function of log10(brake-wear concentration + 1) in the Experiment 3 brake-wear size comparison. Predictions are coloured by brake-wear size class, raw observations are overlaid, 95% confidence intervals are shown, and panels are shown by day from start."
} else {
  "Model-predicted motile fraction as a function of log10(brake-wear concentration + 1) in the Experiment 3 brake-wear size comparison. Predictions are coloured by brake-wear size class, raw observations are overlaid, and 95% confidence intervals are shown. Experiment 3 was sampled on a single day only, so no day facet is shown."
}

if (exists("write_figure_caption_md", mode = "function", inherits = TRUE)) {
  write_figure_caption_md(
    figure_name = "Fig4_brake_size_model_predictions",
    caption_text = caption_text,
    subdir = exp3_figure_subdir
  )
} else {
  write_caption_md_local(
    figure_name = "Fig4_brake_size_model_predictions",
    caption_text = caption_text,
    subdir = exp3_figure_subdir
  )
}

cat("Experiment 3 figure processed.\n\n")

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
  lrt_additive_vs_full = lrt_additive_vs_full,
  lrt_partial_vs_full = lrt_partial_vs_full,
  lrt_null_vs_additive = lrt_null_vs_additive,
  preferred_model_name = preferred_model_name
)

if (exists("write_model_summary_rds", mode = "function", inherits = TRUE)) {
  write_model_summary_rds(
    object = selection_bundle,
    file_stem = "exp3_brake_size_model_tests_and_selection",
    subdir = exp3_model_subdir,
    quiet = FALSE
  )
  
  write_model_summary_rds(
    object = pred_out,
    file_stem = "exp3_brake_size_preferred_model_predictions",
    subdir = exp3_model_subdir,
    quiet = FALSE
  )
} else {
  save_rds_local(
    object = selection_bundle,
    file_path = file.path(dir_model_output, "exp3_brake_size_model_tests_and_selection.rds")
  )
  
  save_rds_local(
    object = pred_out,
    file_path = file.path(dir_model_output, "exp3_brake_size_preferred_model_predictions.rds")
  )
}

# ---------------------------------------------------------
# 20. Optional emmeans summaries
# ---------------------------------------------------------

emmeans_size <- NULL
emmeans_day <- NULL

# Only calculate simple marginal means when the preferred
# model does not rely on a size-by-concentration interaction
# as the primary interpretive result.
preferred_has_size_conc_interaction <- preferred_model_name == "full_primary"

if (exists("get_emmeans_table", mode = "function", inherits = TRUE) &&
    !preferred_has_size_conc_interaction) {
  
  emmeans_size <- tryCatch(
    get_emmeans_table(
      model = preferred_model,
      specs = ~ particle_size_class,
      type = "response"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(emmeans_size) && exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
    write_model_summary_csv(
      data = emmeans_size,
      file_stem = "exp3_brake_size_emmeans_by_size_class",
      subdir = exp3_model_subdir,
      quiet = FALSE
    )
  }
  
  if (has_day_variation) {
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
        file_stem = "exp3_brake_size_emmeans_by_day",
        subdir = exp3_model_subdir,
        quiet = FALSE
      )
    }
  }
}

# ---------------------------------------------------------
# 21. Write log
# ---------------------------------------------------------

sink(file_model_log)
cat("SPM Analysis - 08_models_exp3_brake_size log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")

cat("Input dataset:\n")
cat(input_source_used, "\n\n")

cat("Experiment 3 subset dimensions:\n")
cat("Rows:", nrow(exp3_dat), "\n")
cat("Columns:", ncol(exp3_dat), "\n\n")

cat("Observed size levels:\n")
cat(paste(observed_size_levels, collapse = ", "), "\n\n")

cat("Observed days:\n")
if (length(observed_days) == 0) {
  cat("No non-missing day values available.\n\n")
} else {
  cat(paste(observed_days, collapse = ", "), "\n\n")
}

cat("Observed brake-wear concentration range (ug/L):\n")
cat(
  paste0(
    round(min(exp3_dat$bw_ug_l, na.rm = TRUE), 3),
    " to ",
    round(max(exp3_dat$bw_ug_l, na.rm = TRUE), 3)
  ),
  "\n\n"
)

cat("Detected design features:\n")
cat("has_day_variation:", has_day_variation, "\n")
cat("has_culture_variation:", has_culture_variation, "\n")
cat("n_culture_levels:", n_culture_levels, "\n\n")

cat("Candidate formulas used:\n")
for (nm in names(candidate_formula_strings)) {
  cat("-", nm, ":", candidate_formula_strings[[nm]], "\n")
}
cat("\n")

cat("Experiment 3 dataset summary:\n")
print(exp3_dataset_summary)
cat("\n")

cat("Model comparison:\n")
print(model_comparison)
cat("\n")

cat("Preferred model name:\n")
cat(preferred_model_name, "\n\n")

cat("Preferred fixed effects:\n")
print(preferred_fixed_effects)
cat("\n")

cat("Likelihood ratio test: additive vs full\n")
print(lrt_additive_vs_full)
cat("\n")

cat("Likelihood ratio test: partial vs full\n")
print(lrt_partial_vs_full)
cat("\n")

cat("Likelihood ratio test: null vs additive\n")
print(lrt_null_vs_additive)
cat("\n")

if (!is.null(emmeans_size)) {
  cat("Estimated marginal means by size class:\n")
  print(emmeans_size)
  cat("\n")
}

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
# 22. Console summary
# ---------------------------------------------------------

cat("Model comparison:\n")
print(model_comparison)
cat("\n")

cat("Preferred model:\n")
cat(preferred_model_name, "\n\n")

cat("Preferred fixed effects:\n")
print(preferred_fixed_effects)
cat("\n")

cat("08 Experiment 3 modelling outputs written successfully.\n\n")

cat("========================================================\n")
cat("SCRIPT 08 COMPLETE: MODELS EXP3 BRAKE SIZE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")