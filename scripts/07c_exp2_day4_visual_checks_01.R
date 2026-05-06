# =========================================================
# Script title: 07c_exp2_day4_visual_checks_01.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4 / Article 4
# Purpose: Generate exploratory and supplementary visual
#          checks for the Day 4 Experiment 2 defined-particle
#          model outputs. This script is intended to support
#          figure selection and biological interpretation,
#          not to refit models.
#
# Inputs:
# - data_derived/day4/ms4_day4_analysis_derived.rds
#   or data_derived/day4/ms4_day4_analysis_derived.csv
# - outputs/tables/day4/07b_exp2_defined_particles_day4_predictions.csv
# - outputs/tables/day4/07b_exp2_defined_particles_day4_relative_effects.csv
# - outputs/tables/day4/07b_exp2_defined_particles_day4_model_comparison.csv
#
# Outputs:
# - outputs/figures/day4/Fig5_exp2_defined_particles_day4_faceted_candidate.{pdf,png,tiff}
# - outputs/figures/day4/FigS_exp2_defined_particles_day4_raw_summary.{pdf,png,tiff}
# - outputs/tables/day4/07c_exp2_day4_visual_check_raw_summary.csv
# - outputs/tables/day4/07c_exp2_day4_visual_check_manifest.csv
# - outputs/logs/day4/07c_exp2_day4_visual_checks_log_*.txt
#
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 04b_derive_day4_analysis_dataset.R before this script.
# - Run 07b_models_exp2_defined_particles_day4_01.R before this script.
# - This script does not change model objects or statistical outputs.
# - The preferred Experiment 2 model selected by Script 07b was
#   particle_only_betabinomial, meaning model predictions are
#   expected to be constant across NTU within each particle type.
# - Main visual aim:
#   Show whether the saved model predictions align with the raw
#   Day 4 well-level observations for sand, kaolinite, and peat.
# - Figure 8 reference condition for Experiment 2:
#   Sand at 25 NTU.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 07c: EXP2 DAY 4 VISUAL CHECKS\n")
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
  "tibble"
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
  "helpers_tables.R"
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

file_day4_rds <- file.path(
  dir_data_derived,
  "day4",
  "ms4_day4_analysis_derived.rds"
)

file_day4_csv <- file.path(
  dir_data_derived,
  "day4",
  "ms4_day4_analysis_derived.csv"
)

dir_day4_tables <- file.path(project_root, "outputs", "tables", "day4")
dir_day4_logs <- file.path(project_root, "outputs", "logs", "day4")
dir_day4_figures <- file.path(project_root, "outputs", "figures", "day4")

dir.create(dir_day4_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_figures, recursive = TRUE, showWarnings = FALSE)

file_predictions <- file.path(
  dir_day4_tables,
  "07b_exp2_defined_particles_day4_predictions.csv"
)

file_relative_effects <- file.path(
  dir_day4_tables,
  "07b_exp2_defined_particles_day4_relative_effects.csv"
)

file_model_comparison <- file.path(
  dir_day4_tables,
  "07b_exp2_defined_particles_day4_model_comparison.csv"
)

file_raw_summary <- file.path(
  dir_day4_tables,
  "07c_exp2_day4_visual_check_raw_summary.csv"
)

file_manifest <- file.path(
  dir_day4_tables,
  "07c_exp2_day4_visual_check_manifest.csv"
)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_log <- file.path(
  dir_day4_logs,
  paste0("07c_exp2_day4_visual_checks_log_", timestamp_now, ".txt")
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

safe_se <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
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

save_plot_multi_local <- function(plot,
                                  file_stem,
                                  width = 8.5,
                                  height = 5.8,
                                  dpi = 600) {
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

if (file.exists(file_day4_rds)) {
  day4_dat <- readRDS(file_day4_rds)
  day4_source_used <- file_day4_rds
} else if (file.exists(file_day4_csv)) {
  day4_dat <- readr::read_csv(file_day4_csv, show_col_types = FALSE)
  day4_source_used <- file_day4_csv
} else {
  stop(
    paste0(
      "No Day 4 derived dataset found.\nExpected one of:\n- ",
      file_day4_rds,
      "\n- ",
      file_day4_csv,
      "\nPlease run 04b_derive_day4_analysis_dataset.R first."
    ),
    call. = FALSE
  )
}

cat("Loaded Day 4 dataset from:\n")
cat(day4_source_used, "\n")
cat("Rows:", nrow(day4_dat), "\n")
cat("Columns:", ncol(day4_dat), "\n\n")

# ---------------------------------------------------------
# 7. Load saved Experiment 2 model outputs
# ---------------------------------------------------------

if (!file.exists(file_predictions)) {
  stop(
    paste0(
      "Prediction file not found:\n",
      file_predictions,
      "\nPlease run 07b_models_exp2_defined_particles_day4_01.R first."
    ),
    call. = FALSE
  )
}

predictions <- readr::read_csv(
  file_predictions,
  show_col_types = FALSE
)

cat("Loaded Experiment 2 predictions from:\n")
cat(file_predictions, "\n")
cat("Rows:", nrow(predictions), "\n\n")

if (file.exists(file_relative_effects)) {
  relative_effects <- readr::read_csv(
    file_relative_effects,
    show_col_types = FALSE
  )
  cat("Loaded Experiment 2 relative effects from:\n")
  cat(file_relative_effects, "\n")
  cat("Rows:", nrow(relative_effects), "\n\n")
} else {
  relative_effects <- NULL
  cat("Relative-effects file not found. Continuing without it.\n\n")
}

if (file.exists(file_model_comparison)) {
  model_comparison <- readr::read_csv(
    file_model_comparison,
    show_col_types = FALSE
  )
  cat("Loaded Experiment 2 model comparison from:\n")
  cat(file_model_comparison, "\n")
  cat("Rows:", nrow(model_comparison), "\n\n")
} else {
  model_comparison <- NULL
  cat("Model-comparison file not found. Continuing without it.\n\n")
}

# ---------------------------------------------------------
# 8. Check required columns
# ---------------------------------------------------------

required_day4_columns <- c(
  "experiment_num",
  "days_from_start",
  "mobile_cell_count",
  "stationary_cell_count",
  "total_cells",
  "motility_ratio",
  "ntu",
  "culture",
  "well"
)

missing_day4_columns <- required_day4_columns[
  !required_day4_columns %in% names(day4_dat)
]

if (length(missing_day4_columns) > 0) {
  stop(
    paste0(
      "Day 4 dataset is missing required column(s):\n- ",
      paste(missing_day4_columns, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

required_prediction_columns <- c(
  "particle_type_model",
  "particle_type_plot",
  "ntu",
  "fit_response",
  "conf_low_response",
  "conf_high_response",
  "model"
)

missing_prediction_columns <- required_prediction_columns[
  !required_prediction_columns %in% names(predictions)
]

if (length(missing_prediction_columns) > 0) {
  stop(
    paste0(
      "Prediction file is missing required column(s):\n- ",
      paste(missing_prediction_columns, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

cat("Required columns verified.\n\n")

# ---------------------------------------------------------
# 9. Recreate Experiment 2 Day 4 raw plotting dataset
# ---------------------------------------------------------

exp2_raw <- day4_dat |>
  dplyr::mutate(
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    ntu = suppressWarnings(as.numeric(ntu)),
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
    experiment_num == 10.2,
    days_from_start == 4,
    !is.na(ntu),
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    total_cells > 0
  )

if (nrow(exp2_raw) == 0) {
  stop("No valid Day 4 Experiment 2 raw rows found.", call. = FALSE)
}

# Resolve particle identity from available columns.
if ("particle_type_plot" %in% names(exp2_raw)) {
  exp2_raw$particle_type_source <- exp2_raw$particle_type_plot
} else if ("particle_type" %in% names(exp2_raw)) {
  exp2_raw$particle_type_source <- exp2_raw$particle_type
} else if ("toxin_exposure" %in% names(exp2_raw)) {
  exp2_raw$particle_type_source <- exp2_raw$toxin_exposure
} else {
  stop(
    "No particle identity column found. Expected particle_type_plot, particle_type, or toxin_exposure.",
    call. = FALSE
  )
}

exp2_raw <- exp2_raw |>
  dplyr::mutate(
    particle_type_model = standardise_upper_trim(particle_type_source),
    particle_type_model = dplyr::case_when(
      particle_type_model %in% c("SAND") ~ "SAND",
      particle_type_model %in% c("KAOLINITE", "KAOLIN") ~ "KAOLINITE",
      particle_type_model %in% c("PEAT") ~ "PEAT",
      particle_type_model %in% c("CONTROL", "NO PARTICLES", "NO_PARTICLES") ~ "CONTROL",
      TRUE ~ particle_type_model
    )
  ) |>
  dplyr::filter(
    particle_type_model %in% c("SAND", "KAOLINITE", "PEAT")
  ) |>
  dplyr::mutate(
    particle_type_model = factor(
      particle_type_model,
      levels = c("SAND", "KAOLINITE", "PEAT")
    ),
    particle_type_plot = factor(
      particle_type_model,
      levels = c("SAND", "KAOLINITE", "PEAT"),
      labels = c("Sand", "Kaolinite", "Peat")
    ),
    ntu_f = factor(
      ntu,
      levels = c(25, 100, 400),
      labels = c("25 NTU", "100 NTU", "400 NTU")
    ),
    culture = as.factor(culture),
    well = as.factor(well)
  )

if (nrow(exp2_raw) == 0) {
  stop("No valid SAND/KAOLINITE/PEAT rows remain for Experiment 2 visual checks.", call. = FALSE)
}

cat("Experiment 2 raw Day 4 plotting dataset created.\n")
cat("Rows:", nrow(exp2_raw), "\n")
cat("Particle types:", paste(levels(droplevels(exp2_raw$particle_type_model)), collapse = ", "), "\n")
cat("NTU values:", paste(sort(unique(exp2_raw$ntu)), collapse = ", "), "\n")
cat("Cultures:", safe_n_distinct(exp2_raw$culture), "\n")
cat("Wells:", safe_n_distinct(exp2_raw$well), "\n\n")

# ---------------------------------------------------------
# 10. Prepare prediction plotting dataset
# ---------------------------------------------------------

plot_pred <- predictions |>
  dplyr::mutate(
    ntu = suppressWarnings(as.numeric(ntu)),
    particle_type_model = standardise_upper_trim(particle_type_model),
    particle_type_model = factor(
      particle_type_model,
      levels = c("SAND", "KAOLINITE", "PEAT")
    ),
    particle_type_plot = as.character(particle_type_plot),
    particle_type_plot = dplyr::case_when(
      particle_type_plot %in% c("SAND", "Sand") ~ "Sand",
      particle_type_plot %in% c("KAOLINITE", "Kaolinite") ~ "Kaolinite",
      particle_type_plot %in% c("PEAT", "Peat") ~ "Peat",
      TRUE ~ particle_type_plot
    ),
    particle_type_plot = factor(
      particle_type_plot,
      levels = c("Sand", "Kaolinite", "Peat")
    ),
    ntu_f = factor(
      ntu,
      levels = c(25, 100, 400),
      labels = c("25 NTU", "100 NTU", "400 NTU")
    )
  )

if (any(is.na(plot_pred$particle_type_plot))) {
  warning(
    "Some prediction rows have missing particle_type_plot after standardisation.",
    call. = FALSE
  )
}

if (any(is.na(plot_pred$ntu_f))) {
  warning(
    "Some prediction rows have NTU values outside 25, 100, 400.",
    call. = FALSE
  )
}

cat("Prediction plotting dataset prepared.\n")
cat("Rows:", nrow(plot_pred), "\n\n")

# ---------------------------------------------------------
# 11. Raw summary for visual checking
# ---------------------------------------------------------

raw_summary <- exp2_raw |>
  dplyr::group_by(
    particle_type_model,
    particle_type_plot,
    ntu,
    ntu_f
  ) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_cultures = safe_n_distinct(culture),
    n_wells = safe_n_distinct(well),
    n_videos = if ("video_file" %in% names(exp2_raw)) safe_n_distinct(video_file) else NA_integer_,
    total_mobile_cells = sum(mobile_cell_count, na.rm = TRUE),
    total_stationary_cells = sum(stationary_cell_count, na.rm = TRUE),
    total_cells = sum(total_cells, na.rm = TRUE),
    mean_motility_ratio = safe_mean(motility_ratio),
    sd_motility_ratio = safe_sd(motility_ratio),
    se_motility_ratio = safe_se(motility_ratio),
    min_motility_ratio = safe_min(motility_ratio),
    max_motility_ratio = safe_max(motility_ratio),
    .groups = "drop"
  ) |>
  dplyr::arrange(particle_type_model, ntu)

readr::write_csv(raw_summary, file_raw_summary)

cat("Raw visual-check summary written to:\n")
cat(file_raw_summary, "\n\n")

print(raw_summary)

# ---------------------------------------------------------
# 12. Report preferred model, if available
# ---------------------------------------------------------

if (!is.null(model_comparison)) {
  preferred_model_name <- model_comparison |>
    dplyr::filter(model_status == "success", !is.na(AIC)) |>
    dplyr::arrange(AIC) |>
    dplyr::slice(1) |>
    dplyr::pull(model)
  
  cat("Preferred model from saved comparison:\n")
  cat(preferred_model_name, "\n\n")
  
  if (preferred_model_name %in% c("particle_only_betabinomial", "particle_only_binomial")) {
    cat("Note: preferred model is particle-only, so model predictions are expected to be flat across NTU within each particle type.\n\n")
  }
} else {
  preferred_model_name <- NA_character_
}

# ---------------------------------------------------------
# 13. Plot A: candidate main Figure 5, faceted by particle type
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

plot_raw <- exp2_raw |>
  dplyr::mutate(
    particle_type_plot = factor(
      particle_type_plot,
      levels = c("Sand", "Kaolinite", "Peat")
    ),
    ntu_f = factor(
      ntu_f,
      levels = c("25 NTU", "100 NTU", "400 NTU")
    )
  )

plot_pred <- plot_pred |>
  dplyr::mutate(
    particle_type_plot = factor(
      particle_type_plot,
      levels = c("Sand", "Kaolinite", "Peat")
    ),
    ntu_f = factor(
      ntu_f,
      levels = c("25 NTU", "100 NTU", "400 NTU")
    )
  )

p_faceted <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data = plot_raw,
    ggplot2::aes(
      x = ntu_f,
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
      x = ntu_f,
      ymin = conf_low_response,
      ymax = conf_high_response
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    data = plot_pred,
    ggplot2::aes(
      x = ntu_f,
      y = fit_response
    ),
    size = 2.8
  ) +
  ggplot2::facet_wrap(~ particle_type_plot, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Experiment 2: Defined particle concentration series",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Turbidity treatment",
    y = "Motile fraction",
    caption = paste(
      "Points show raw well-level observations.",
      "Black points and error bars show saved model predictions with 95% confidence intervals.",
      "Predictions are expected to be flat across NTU if the preferred model is particle-only."
    )
  )

fig_faceted_manifest <- save_plot_multi_local(
  plot = p_faceted,
  file_stem = "Fig5_exp2_defined_particles_day4_faceted_candidate",
  width = 9,
  height = 5.8,
  dpi = 600
)

cat("Faceted candidate figure saved:\n")
print(fig_faceted_manifest)

# ---------------------------------------------------------
# 14. Plot B: raw-data-only summary / sanity check
# ---------------------------------------------------------

p_raw_summary <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data = plot_raw,
    ggplot2::aes(
      x = ntu_f,
      y = motility_ratio
    ),
    width = 0.08,
    height = 0,
    alpha = 0.35,
    size = 1.8
  ) +
  ggplot2::geom_errorbar(
    data = raw_summary,
    ggplot2::aes(
      x = ntu_f,
      ymin = mean_motility_ratio - se_motility_ratio,
      ymax = mean_motility_ratio + se_motility_ratio
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    data = raw_summary,
    ggplot2::aes(
      x = ntu_f,
      y = mean_motility_ratio
    ),
    size = 2.8
  ) +
  ggplot2::facet_wrap(~ particle_type_plot, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Experiment 2: Raw Day 4 motility ratios",
    subtitle = "Raw well-level observations with group means ± SE",
    x = "Turbidity treatment",
    y = "Motility ratio",
    caption = paste(
      "This plot is a visual sanity check only.",
      "It shows raw group means and standard errors, not model-estimated confidence intervals."
    )
  )

fig_raw_manifest <- save_plot_multi_local(
  plot = p_raw_summary,
  file_stem = "FigS_exp2_defined_particles_day4_raw_summary",
  width = 9,
  height = 5.8,
  dpi = 600
)

cat("Raw summary figure saved:\n")
print(fig_raw_manifest)

# ---------------------------------------------------------
# 15. Optional Plot C: compact non-faceted grouped view
# ---------------------------------------------------------

particle_palette_local <- c(
  "Sand" = "#C2B280",
  "Kaolinite" = "#A6CEE3",
  "Peat" = "#8B4513"
)

p_grouped <- ggplot2::ggplot() +
  ggplot2::geom_jitter(
    data = plot_raw,
    ggplot2::aes(
      x = ntu_f,
      y = motility_ratio,
      colour = particle_type_plot
    ),
    width = 0.10,
    height = 0,
    alpha = 0.35,
    size = 1.8
  ) +
  ggplot2::geom_errorbar(
    data = plot_pred,
    ggplot2::aes(
      x = ntu_f,
      ymin = conf_low_response,
      ymax = conf_high_response,
      colour = particle_type_plot
    ),
    width = 0.15,
    linewidth = 0.6,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::geom_point(
    data = plot_pred,
    ggplot2::aes(
      x = ntu_f,
      y = fit_response,
      colour = particle_type_plot
    ),
    size = 2.8,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::scale_colour_manual(values = particle_palette_local, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Experiment 2: Defined particle concentration series",
    subtitle = "Compact grouped view of saved Day 4 model predictions",
    x = "Turbidity treatment",
    y = "Motile fraction",
    colour = "Particle type",
    caption = paste(
      "Points show raw observations.",
      "Model predictions and 95% confidence intervals are overlaid by particle type."
    )
  )

fig_grouped_manifest <- save_plot_multi_local(
  plot = p_grouped,
  file_stem = "Fig5_exp2_defined_particles_day4_grouped_candidate",
  width = 8,
  height = 5.8,
  dpi = 600
)

cat("Grouped candidate figure saved:\n")
print(fig_grouped_manifest)

# ---------------------------------------------------------
# 16. Write figure manifest
# ---------------------------------------------------------

visual_manifest <- dplyr::bind_rows(
  fig_faceted_manifest |>
    dplyr::mutate(
      figure_role = "candidate_main_figure",
      description = "Faceted by particle type; raw points plus saved model predictions and 95% confidence intervals."
    ),
  fig_raw_manifest |>
    dplyr::mutate(
      figure_role = "supplementary_or_sanity_check",
      description = "Raw well-level data with group means ± SE; no model predictions."
    ),
  fig_grouped_manifest |>
    dplyr::mutate(
      figure_role = "alternative_compact_view",
      description = "Grouped particle-colour view; raw points plus saved model predictions and 95% confidence intervals."
    )
)

readr::write_csv(visual_manifest, file_manifest)

cat("Visual-check manifest written to:\n")
cat(file_manifest, "\n\n")

print(visual_manifest)

# ---------------------------------------------------------
# 17. Write log
# ---------------------------------------------------------

sink(file_log)
cat("Experiment 2 Day 4 visual checks log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Day 4 data source:\n")
cat(day4_source_used, "\n\n")

cat("Prediction file:\n")
cat(file_predictions, "\n\n")

cat("Relative-effects file:\n")
cat(ifelse(file.exists(file_relative_effects), file_relative_effects, "Not found"), "\n\n")

cat("Model-comparison file:\n")
cat(ifelse(file.exists(file_model_comparison), file_model_comparison, "Not found"), "\n\n")

cat("Rows in Experiment 2 raw visual-check dataset:", nrow(exp2_raw), "\n")
cat("Particle types:\n")
print(table(exp2_raw$particle_type_model, useNA = "ifany"))

cat("\nNTU values:\n")
print(sort(unique(exp2_raw$ntu)))

cat("\nCultures:", safe_n_distinct(exp2_raw$culture), "\n")
cat("Wells:", safe_n_distinct(exp2_raw$well), "\n\n")

cat("Preferred model from saved model comparison:\n")
cat(preferred_model_name, "\n\n")

cat("Raw visual-check summary:\n")
print(raw_summary)

cat("\nFigure manifest:\n")
print(visual_manifest)

sink()

cat("Log written to:\n")
cat(file_log, "\n\n")

cat("SCRIPT 07c COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
