# =========================================================
# Script title: 07_models_exp2_defined_particles_v3.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Fit and evaluate Experiment 2 defined-particle
#          concentration-series models for motile fraction
#          using primary cell counts, generate prediction
#          outputs, compare candidate models, and save
#          experiment-specific summaries and figures to disk.
# Inputs:
# - data_derived/ms4_analysis_derived.csv or .rds
# Outputs:
# - outputs/models/exp2_defined_particles/*.rds
# - outputs/models/exp2_defined_particles/*_fixed_effects.csv
# - outputs/models/exp2_defined_particles/*_metadata.csv
# - outputs/models/exp2_defined_particles/*_diagnostics_summary.csv
# - outputs/models/exp2_defined_particles/*_model_comparison.csv
# - outputs/models/exp2_defined_particles/*_predictions.csv
# - outputs/tables/07_exp2_defined_particles_dataset_summary.csv
# - outputs/tables/07_exp2_defined_particles_model_comparison.csv
# - outputs/tables/07_exp2_defined_particles_model_summary_table.csv
# - outputs/figures/models_exp2/Fig3_defined_particles_model_predictions_*.{pdf,png,tiff}
# - outputs/logs/07_models_exp2_defined_particles_log_*.txt
# Date created: 26 March 2026
# Last updated: 26 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 04_derive_variables.R before this script.
# - Primary fitted response:
#   cbind(mobile_cell_count, stationary_cell_count)
# - Manuscript 4 block must be restricted to:
#   experiment_num %in% c(8.2, 9.2, 10.2, 11.2)
# - Focal modelling subset for this script:
#   experiment_num == 10.2
# - Core candidate model from brief:
#   cbind(mobile_cell_count, stationary_cell_count) ~
#   log10(ntu + 1) * particle_type + days_from_start + (1 | culture)
# - Experiment 2 particle classes:
#   Sand, Kaolinite, Peat
# - Controls may exist within the same experimental block,
#   but treatment interpretation is only valid inside the
#   focal experiment subset.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 07: MODELS EXP2 DEFINED PARTICLES\n")
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
# 3. Define input and output paths
# ---------------------------------------------------------

file_input_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")

exp2_model_subdir <- "exp2_defined_particles"
exp2_figure_subdir <- "models_exp2"

file_dataset_summary <- file.path(dir_tables, "07_exp2_defined_particles_dataset_summary.csv")
file_model_comparison_table <- file.path(dir_tables, "07_exp2_defined_particles_model_comparison.csv")
file_model_summary_table <- file.path(dir_tables, "07_exp2_defined_particles_model_summary_table.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_model_log <- file.path(
  dir_logs,
  paste0("07_models_exp2_defined_particles_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 4. Import derived dataset
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
# 5. Check required columns
# ---------------------------------------------------------

required_columns <- c(
  "experiment",
  "experiment_num",
  "experiment_type",
  "mobile_cell_count",
  "stationary_cell_count",
  "ntu",
  "days_from_start",
  "culture"
)

missing_required_columns <- required_columns[!required_columns %in% names(dat)]

if (length(missing_required_columns) > 0) {
  stop(
    paste0(
      "The following required Experiment 2 column(s) are missing:\n- ",
      paste(missing_required_columns, collapse = "\n- "),
      "\nPlease review upstream scripts."
    ),
    call. = FALSE
  )
}

cat("Required Experiment 2 columns verified.\n\n")

# ---------------------------------------------------------
# 6. Restrict to Manuscript 4 block, then focal experiment
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

exp2_dat <- dat |>
  dplyr::filter(experiment_num == 10.2) |>
  dplyr::filter(
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    !is.na(ntu),
    !is.na(days_from_start),
    !is.na(culture)
  ) |>
  dplyr::mutate(
    culture = as.factor(culture),
    days_from_start = as.integer(days_from_start),
    ntu = as.numeric(ntu)
  )

if (nrow(exp2_dat) == 0) {
  stop("No Experiment 2 rows available after filtering.", call. = FALSE)
}

cat("Experiment 2 focal subset created.\n")
cat("Rows:", nrow(exp2_dat), "\n")
cat("Cultures:", dplyr::n_distinct(exp2_dat$culture, na.rm = TRUE), "\n")
cat("Days:", paste(sort(unique(exp2_dat$days_from_start)), collapse = ", "), "\n")
cat("NTU values observed:", paste(sort(unique(exp2_dat$ntu)), collapse = ", "), "\n\n")

# ---------------------------------------------------------
# 7. Standardise particle identity and derive support vars
# ---------------------------------------------------------

standardise_upper_trim <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
}

if ("particle_type_plot" %in% names(exp2_dat)) {
  exp2_dat$particle_type_model <- standardise_upper_trim(exp2_dat$particle_type_plot)
} else if ("particle_type" %in% names(exp2_dat)) {
  exp2_dat$particle_type_model <- standardise_upper_trim(exp2_dat$particle_type)
} else if ("toxin_exposure" %in% names(exp2_dat)) {
  exp2_dat$particle_type_model <- standardise_upper_trim(exp2_dat$toxin_exposure)
} else {
  exp2_dat$particle_type_model <- NA_character_
}

exp2_dat <- exp2_dat |>
  dplyr::mutate(
    particle_type_model = dplyr::case_when(
      particle_type_model %in% c("SAND") ~ "SAND",
      particle_type_model %in% c("KAOLINITE") ~ "KAOLINITE",
      particle_type_model %in% c("PEAT") ~ "PEAT",
      particle_type_model %in% c("CONTROL") ~ "CONTROL",
      TRUE ~ particle_type_model
    ),
    total_cells = mobile_cell_count + stationary_cell_count,
    motility_ratio = dplyr::if_else(
      total_cells > 0,
      mobile_cell_count / total_cells,
      NA_real_
    ),
    log10_ntu_plus1 = log10(ntu + 1),
    days_from_start_f = factor(
      days_from_start,
      levels = sort(unique(days_from_start))
    ),
    experiment_plot_label = "Experiment 2: Defined particles"
  )

# Restrict to focal particle classes for primary interpretation.
# Retain CONTROL rows only if present in exp2 itself.
exp2_dat <- exp2_dat |>
  dplyr::filter(!is.na(particle_type_model)) |>
  dplyr::mutate(
    particle_type_model = factor(
      particle_type_model,
      levels = c("CONTROL", "SAND", "KAOLINITE", "PEAT")
    )
  )

observed_particle_levels <- levels(droplevels(exp2_dat$particle_type_model))

cat("Experiment 2 modelling support variables derived.\n")
cat("Particle types observed:", paste(observed_particle_levels, collapse = ", "), "\n\n")

# ---------------------------------------------------------
# 8. Build dataset summary table
# ---------------------------------------------------------

has_well <- "well" %in% names(exp2_dat)
has_video_file <- "video_file" %in% names(exp2_dat)

exp2_dataset_summary <- exp2_dat |>
  dplyr::group_by(days_from_start, particle_type_model) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
    n_wells = if (has_well) dplyr::n_distinct(well, na.rm = TRUE) else NA_integer_,
    n_videos = if (has_video_file) dplyr::n_distinct(video_file, na.rm = TRUE) else NA_integer_,
    min_ntu = min(ntu, na.rm = TRUE),
    max_ntu = max(ntu, na.rm = TRUE),
    mean_total_cells = mean(total_cells, na.rm = TRUE),
    mean_motility_ratio = mean(motility_ratio, na.rm = TRUE),
    sd_motility_ratio = stats::sd(motility_ratio, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(days_from_start, particle_type_model)

readr::write_csv(exp2_dataset_summary, file_dataset_summary)

cat("Experiment 2 dataset summary written to:\n")
cat(file_dataset_summary, "\n\n")

# ---------------------------------------------------------
# 9. Define candidate model formulas
# ---------------------------------------------------------

formula_full_primary <- cbind(mobile_cell_count, stationary_cell_count) ~
  log10_ntu_plus1 * particle_type_model * days_from_start + (1 | culture)

formula_additive_primary <- cbind(mobile_cell_count, stationary_cell_count) ~
  log10_ntu_plus1 + particle_type_model + days_from_start + (1 | culture)

formula_partial_primary <- cbind(mobile_cell_count, stationary_cell_count) ~
  log10_ntu_plus1 * particle_type_model + days_from_start + (1 | culture)

formula_null <- cbind(mobile_cell_count, stationary_cell_count) ~
  1 + (1 | culture)

cat("Candidate model formulas defined.\n\n")

# ---------------------------------------------------------
# 10. Fit candidate models
# ---------------------------------------------------------

if (!exists("fit_check_save_model", mode = "function", inherits = TRUE)) {
  stop("helpers_model_checks.R was not loaded correctly.", call. = FALSE)
}

fit_full_primary <- fit_check_save_model(
  formula = formula_full_primary,
  data = exp2_dat,
  model_name = "exp2_defined_particles_full_primary",
  model_subdir = exp2_model_subdir,
  save_terms_csv = TRUE,
  save_meta_csv = TRUE,
  save_diagnostics = TRUE,
  run_dharma = TRUE,
  run_performance = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  quiet = FALSE
)

fit_additive_primary <- fit_check_save_model(
  formula = formula_additive_primary,
  data = exp2_dat,
  model_name = "exp2_defined_particles_additive_primary",
  model_subdir = exp2_model_subdir,
  save_terms_csv = TRUE,
  save_meta_csv = TRUE,
  save_diagnostics = TRUE,
  run_dharma = TRUE,
  run_performance = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  quiet = FALSE
)

fit_partial_primary <- fit_check_save_model(
  formula = formula_partial_primary,
  data = exp2_dat,
  model_name = "exp2_defined_particles_partial_primary",
  model_subdir = exp2_model_subdir,
  save_terms_csv = TRUE,
  save_meta_csv = TRUE,
  save_diagnostics = TRUE,
  run_dharma = TRUE,
  run_performance = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  quiet = FALSE
)

fit_null <- fit_check_save_model(
  formula = formula_null,
  data = exp2_dat,
  model_name = "exp2_defined_particles_null",
  model_subdir = exp2_model_subdir,
  save_terms_csv = TRUE,
  save_meta_csv = TRUE,
  save_diagnostics = TRUE,
  run_dharma = TRUE,
  run_performance = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  quiet = FALSE
)

cat("Candidate models fitted.\n\n")

# check
if (!all(c("fit_full_primary", "fit_additive_primary", "fit_partial_primary", "fit_null") %in% ls())) {
  stop("One or more Experiment 2 model fits failed. Check errors above.", call. = FALSE)
}

required_models <- c(
  "fit_full_primary",
  "fit_additive_primary",
  "fit_partial_primary",
  "fit_null"
)

missing_models <- required_models[!sapply(required_models, exists)]

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
# 11. Compare candidate models
# ---------------------------------------------------------

if (!exists("compare_models_aic", mode = "function", inherits = TRUE)) {
  stop("compare_models_aic() is not available.", call. = FALSE)
}

model_comparison <- compare_models_aic(
  full_primary = fit_full_primary$model,
  additive_primary = fit_additive_primary$model,
  partial_primary = fit_partial_primary$model,
  null = fit_null$model
) |>
  dplyr::mutate(
    experiment = "Experiment 2",
    experiment_num = 10.2
  )

readr::write_csv(model_comparison, file_model_comparison_table)

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = model_comparison,
    file_stem = "exp2_defined_particles_model_comparison",
    subdir = exp2_model_subdir,
    quiet = FALSE
  )
}

cat("Model comparison written to:\n")
cat(file_model_comparison_table, "\n\n")

# ---------------------------------------------------------
# 12. Likelihood ratio comparisons
# ---------------------------------------------------------

lrt_additive_vs_full <- tryCatch(
  stats::anova(fit_additive_primary$model, fit_full_primary$model),
  error = function(e) e
)

lrt_full_vs_day_interaction <- tryCatch(
  stats::anova(fit_additive_primary$model, fit_partial_primary$model),
  error = function(e) e
)

lrt_null_vs_additive <- tryCatch(
  stats::anova(fit_null$model, fit_additive_primary$model),
  error = function(e) e
)

# ---------------------------------------------------------
# 13. Select preferred model
# ---------------------------------------------------------

preferred_model_name <- model_comparison$model[[1]]

preferred_model <- switch(
  preferred_model_name,
  full_primary = fit_full_primary$model,
  additive_primary = fit_additive_primary$model,
  interaction_day_sensitivity = fit_partial_primary$model,
  null = fit_null$model,
  fit_full_primary$model
)

cat("Preferred model selected by lowest AIC:\n")
cat(preferred_model_name, "\n\n")

# ---------------------------------------------------------
# 14. Build prediction grid
# ---------------------------------------------------------

if (!exists("make_prediction_grid", mode = "function", inherits = TRUE) ||
    !exists("predict_glmmtmb_response", mode = "function", inherits = TRUE)) {
  stop("Prediction helpers are not available.", call. = FALSE)
}

pred_grid <- make_prediction_grid(
  data = exp2_dat,
  focal_terms = c("log10_ntu_plus1", "particle_type_model", "days_from_start"),
  at = list(
    particle_type_model = levels(droplevels(exp2_dat$particle_type_model)),
    days_from_start = sort(unique(exp2_dat$days_from_start))
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
    prediction_scale = "log10_ntu_plus1"
  )

if (exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
  write_model_summary_csv(
    data = pred_out,
    file_stem = "exp2_defined_particles_preferred_model_predictions",
    subdir = exp2_model_subdir,
    quiet = FALSE
  )
}

cat("Prediction outputs generated and saved.\n\n")

# ---------------------------------------------------------
# 15. Build fixed-effect summary for preferred model
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
    experiment = "Experiment 2",
    experiment_num = 10.2
  )

readr::write_csv(preferred_fixed_effects, file_model_summary_table)

cat("Preferred model fixed-effect summary written to:\n")
cat(file_model_summary_table, "\n\n")

# ---------------------------------------------------------
# 16. Build Experiment 2 figure
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

if (exists("particle_palette", mode = "function", inherits = TRUE)) {
  exp2_particle_palette <- particle_palette()
} else if (exists("particle_palette", inherits = TRUE)) {
  exp2_particle_palette <- get("particle_palette", inherits = TRUE)
} else {
  exp2_particle_palette <- c(
    "CONTROL"   = "#4D4D4D",
    "SAND"      = "#C2B280",
    "KAOLINITE" = "#A6CEE3",
    "PEAT"      = "#8B4513"
  )
}

plot_raw <- exp2_dat |>
  dplyr::mutate(
    day_label = paste0("Day ", days_from_start)
  )

p_exp2 <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = plot_raw,
    ggplot2::aes(
      x = log10_ntu_plus1,
      y = motility_ratio,
      colour = particle_type_model
    ),
    alpha = 0.30,
    size = 1.8
  ) +
  ggplot2::geom_ribbon(
    data = pred_out,
    ggplot2::aes(
      x = log10_ntu_plus1,
      ymin = conf_low_response,
      ymax = conf_high_response,
      fill = particle_type_model
    ),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = pred_out,
    ggplot2::aes(
      x = log10_ntu_plus1,
      y = fit_response,
      colour = particle_type_model
    ),
    linewidth = 0.8
  ) +
  ggplot2::facet_wrap(~ day_label) +
  ggplot2::scale_colour_manual(values = exp2_particle_palette, drop = FALSE) +
  ggplot2::scale_fill_manual(values = exp2_particle_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Defined particle concentration series",
    subtitle = "Model-predicted motile fraction versus log10(NTU + 1), by particle type and day",
    x = "log10(NTU + 1)",
    y = "Motile fraction",
    colour = "Particle type",
    fill = "Particle type",
    caption = "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals. Panels are shown by day from start."
  )

# ---------------------------------------------------------
# 17. Save Experiment 2 figure
# ---------------------------------------------------------

saved_figures <- character(0)

if (exists("save_figure_both_widths", mode = "function", inherits = TRUE)) {
  saved_figures <- save_figure_both_widths(
    plot = p_exp2,
    figure_name = "Fig3_defined_particles_model_predictions",
    subdir = exp2_figure_subdir,
    height = "standard",
    quiet = TRUE
  )
}

if (exists("write_figure_caption_md", mode = "function", inherits = TRUE)) {
  write_figure_caption_md(
    figure_name = "Fig3_defined_particles_model_predictions",
    caption_text = "Model-predicted motile fraction as a function of log10(NTU + 1) in the defined-particle concentration experiment. Predictions are coloured by particle type, raw observations are overlaid, and 95% confidence intervals are shown. Panels are shown by day from start.",
    subdir = exp2_figure_subdir
  )
}

cat("Experiment 2 figure processed.\n\n")

# ---------------------------------------------------------
# 18. Save extra model objects and summaries
# ---------------------------------------------------------

if (exists("write_model_summary_rds", mode = "function", inherits = TRUE)) {
  write_model_summary_rds(
    object = list(
      lrt_additive_vs_full = lrt_additive_vs_full,
      lrt_full_vs_day_interaction = lrt_full_vs_day_interaction,
      lrt_null_vs_additive = lrt_null_vs_additive,
      preferred_model_name = preferred_model_name
    ),
    file_stem = "exp2_defined_particles_model_tests_and_selection",
    subdir = exp2_model_subdir,
    quiet = FALSE
  )

  write_model_summary_rds(
    object = pred_out,
    file_stem = "exp2_defined_particles_preferred_model_predictions",
    subdir = exp2_model_subdir,
    quiet = FALSE
  )
}

# ---------------------------------------------------------
# 19. Optional emmeans summaries
# ---------------------------------------------------------

emmeans_particle <- NULL
emmeans_day <- NULL

if (exists("get_emmeans_table", mode = "function", inherits = TRUE)) {
  emmeans_particle <- tryCatch(
    get_emmeans_table(
      model = preferred_model,
      specs = ~ particle_type_model,
      type = "response"
    ),
    error = function(e) NULL
  )

  if (!is.null(emmeans_particle) && exists("write_model_summary_csv", mode = "function", inherits = TRUE)) {
    write_model_summary_csv(
      data = emmeans_particle,
      file_stem = "exp2_defined_particles_emmeans_by_particle_type",
      subdir = exp2_model_subdir,
      quiet = FALSE
    )
  }

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
      file_stem = "exp2_defined_particles_emmeans_by_day",
      subdir = exp2_model_subdir,
      quiet = FALSE
    )
  }
}

# ---------------------------------------------------------
# 20. Write log
# ---------------------------------------------------------

sink(file_model_log)
cat("SPM Analysis - 07_models_exp2_defined_particles log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")

cat("Input dataset:\n")
cat(input_source_used, "\n\n")

cat("Experiment 2 subset dimensions:\n")
cat("Rows:", nrow(exp2_dat), "\n")
cat("Columns:", ncol(exp2_dat), "\n\n")

cat("Observed particle levels:\n")
cat(paste(observed_particle_levels, collapse = ", "), "\n\n")

cat("Experiment 2 dataset summary:\n")
print(exp2_dataset_summary)
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

cat("Likelihood ratio test: additive vs day-interaction sensitivity\n")
print(lrt_full_vs_day_interaction)
cat("\n")

cat("Likelihood ratio test: null vs additive\n")
print(lrt_null_vs_additive)
cat("\n")

if (!is.null(emmeans_particle)) {
  cat("Estimated marginal means by particle type:\n")
  print(emmeans_particle)
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
# 21. Console summary
# ---------------------------------------------------------

cat("Model comparison:\n")
print(model_comparison)
cat("\n")

cat("Preferred model:\n")
cat(preferred_model_name, "\n\n")

cat("Preferred fixed effects:\n")
print(preferred_fixed_effects)
cat("\n")

cat("07 Experiment 2 modelling outputs written successfully.\n\n")

cat("========================================================\n")
cat("SCRIPT 07 COMPLETE: MODELS EXP2 DEFINED PARTICLES\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

