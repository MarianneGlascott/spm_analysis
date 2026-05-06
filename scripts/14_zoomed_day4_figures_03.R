# =========================================================
# Script title: 14_zoomed_day4_figures_01.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4 / Article 4
# Purpose: Generate alternative zoomed versions of the Day 4
#          main figures (Figures 4-8) using saved prediction
#          and relative-effect outputs. These figures are
#          regenerated for visual comparison only and do not
#          overwrite the existing plots.
#
# Inputs:
# - data_derived/day4/ms4_day4_analysis_derived.rds
# - outputs/tables/day4/06b_exp1_light_day4_predictions.csv
# - outputs/tables/day4/07b_exp2_defined_particles_day4_predictions.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_predictions.csv
# - outputs/tables/day4/09b_exp4_field_spm_day4_predictions.csv
# - outputs/tables/day4/13_figure8_day4_relative_effects_combined.csv
#
# Outputs:
# - outputs/figures/day4/Fig4_exp1_light_day4_model_predictions_zoomed.{pdf,png,tiff}
# - outputs/figures/day4/Fig5_exp2_defined_particles_day4_model_predictions_zoomed.{pdf,png,tiff}
# - outputs/figures/day4/Fig6_exp3_brake_size_day4_model_predictions_zoomed.{pdf,png,tiff}
# - outputs/figures/day4/Fig7_exp4_field_spm_day4_model_predictions_zoomed.{pdf,png,tiff}
# - outputs/figures/day4/Fig8_day4_relative_effects_forest_plot_zoomed.{pdf,png,tiff}
# - outputs/tables/day4/14_zoomed_day4_figures_manifest.csv
# - outputs/logs/day4/14_zoomed_day4_figures_log_*.txt
#
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - This script does not fit models or modify predictions.
# - All zooming is done using coord_cartesian() so the
#   underlying data and model outputs remain unchanged.
# - Zoom settings:
#   Figures 4-7: y-axis zoom to 0.20-0.55
#   Figure 8:    x-axis zoom to -45 to 50
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 14: ZOOMED DAY 4 FIGURES\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check setup objects
# ---------------------------------------------------------

required_objects <- c(
  "project_root",
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
  "forcats",
  "tibble",
  "scales"
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
  "helpers_labels.R",
  "helpers_model_checks.R",
  "helpers_save_figures.R",
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

dir_day4_derived <- file.path(project_root, "data_derived", "day4")
dir_day4_tables  <- file.path(project_root, "outputs", "tables", "day4")
dir_day4_logs    <- file.path(project_root, "outputs", "logs", "day4")
dir_day4_figures <- file.path(project_root, "outputs", "figures", "day4")

dir.create(dir_day4_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_tables, recursive = TRUE, showWarnings = FALSE)

file_day4_rds <- file.path(dir_day4_derived, "ms4_day4_analysis_derived.rds")
file_day4_csv <- file.path(dir_day4_derived, "ms4_day4_analysis_derived.csv")

file_pred_exp1 <- file.path(dir_day4_tables, "06b_exp1_light_day4_predictions.csv")
file_pred_exp2 <- file.path(dir_day4_tables, "07b_exp2_defined_particles_day4_predictions.csv")
file_pred_exp3 <- file.path(dir_day4_tables, "08b_exp3_brake_size_day4_predictions.csv")
file_pred_exp4 <- file.path(dir_day4_tables, "09b_exp4_field_spm_day4_predictions.csv")
file_fig8_combined <- file.path(dir_day4_tables, "13_figure8_day4_relative_effects_combined.csv")

file_manifest <- file.path(dir_day4_tables, "14_zoomed_day4_figures_manifest.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_log <- file.path(
  dir_day4_logs,
  paste0("14_zoomed_day4_figures_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 5. Settings
# ---------------------------------------------------------

zoom_y_limits <- c(0.20, 0.55)
zoom_x_limits_fig8 <- c(-45, 50)

fig_width_std  <- 8.0
fig_height_std <- 6.0
fig_width_f8   <- 9.5
fig_height_f8  <- 8.5
fig_dpi        <- 600

cat("Zoom settings:\n")
cat("Figures 4-7 y-axis:", paste(zoom_y_limits, collapse = " to "), "\n")
cat("Figure 8 x-axis:", paste(zoom_x_limits_fig8, collapse = " to "), "\n\n")

# ---------------------------------------------------------
# 6. Helper functions
# ---------------------------------------------------------

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      paste0("Required input file not found for ", label, ":\n", path),
      call. = FALSE
    )
  }
  readr::read_csv(path, show_col_types = FALSE)
}

load_day4_dataset <- function(file_rds, file_csv) {
  if (file.exists(file_rds)) {
    obj <- readRDS(file_rds)
    return(obj)
  }
  if (file.exists(file_csv)) {
    obj <- readr::read_csv(file_csv, show_col_types = FALSE)
    return(obj)
  }
  stop(
    paste0(
      "No Day 4 dataset found.\nChecked:\n- ",
      file_rds,
      "\n- ",
      file_csv
    ),
    call. = FALSE
  )
}

save_plot_multi_local <- function(plot,
                                  file_stem,
                                  width = 8,
                                  height = 6,
                                  dpi = 600) {
  file_pdf  <- file.path(dir_day4_figures, paste0(file_stem, ".pdf"))
  file_png  <- file.path(dir_day4_figures, paste0(file_stem, ".png"))
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

check_required_columns <- function(df, required_cols, label) {
  missing_cols <- required_cols[!required_cols %in% names(df)]
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        label, " is missing required column(s):\n- ",
        paste(missing_cols, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

warn_if_zoom_clips <- function(values_low, values_high, zoom_limits, label) {
  min_val <- suppressWarnings(min(values_low, na.rm = TRUE))
  max_val <- suppressWarnings(max(values_high, na.rm = TRUE))
  
  if (is.finite(min_val) && min_val < zoom_limits[1]) {
    warning(
      paste0(
        label, ": lower values extend below zoom limit (",
        round(min_val, 4), " < ", zoom_limits[1], ")."
      ),
      call. = FALSE
    )
  }
  
  if (is.finite(max_val) && max_val > zoom_limits[2]) {
    warning(
      paste0(
        label, ": upper values extend above zoom limit (",
        round(max_val, 4), " > ", zoom_limits[2], ")."
      ),
      call. = FALSE
    )
  }
}

# ---------------------------------------------------------
# 7. Load data
# ---------------------------------------------------------

day4_data <- load_day4_dataset(file_day4_rds, file_day4_csv)

pred_exp1 <- read_required_csv(file_pred_exp1, "Experiment 1 predictions")
pred_exp2 <- read_required_csv(file_pred_exp2, "Experiment 2 predictions")
pred_exp3 <- read_required_csv(file_pred_exp3, "Experiment 3 predictions")
pred_exp4 <- read_required_csv(file_pred_exp4, "Experiment 4 predictions")
fig8_combined <- read_required_csv(file_fig8_combined, "Figure 8 combined effects")

cat("Input files loaded successfully.\n")
cat("Day 4 dataset rows:", nrow(day4_data), "\n")
cat("Exp1 prediction rows:", nrow(pred_exp1), "\n")
cat("Exp2 prediction rows:", nrow(pred_exp2), "\n")
cat("Exp3 prediction rows:", nrow(pred_exp3), "\n")
cat("Exp4 prediction rows:", nrow(pred_exp4), "\n")
cat("Figure 8 combined rows:", nrow(fig8_combined), "\n\n")

# ---------------------------------------------------------
# 8. Prepare plot theme
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

# Fallback colours if helper objects are not found
if (!exists("kelp_raw_point_colour", inherits = TRUE)) {
  kelp_raw_point_colour <- "#b8b8b8"
}
if (!exists("kelp_reference_line_colour", inherits = TRUE)) {
  kelp_reference_line_colour <- "#6c6c6c"
}

# ---------------------------------------------------------
# 9. Figure 4 - Experiment 1 zoomed
# ---------------------------------------------------------

raw_exp1 <- day4_data |>
  dplyr::filter(experiment_num == 9.2)

# Flexible label preparation
if (!"lux_plot_f" %in% names(raw_exp1)) {
  if ("lux_exposure" %in% names(raw_exp1)) {
    raw_exp1 <- raw_exp1 |>
      dplyr::mutate(
        lux_plot_f = factor(
          paste0(lux_exposure, " lux"),
          levels = c("0 lux", "4 lux", "70 lux", "117 lux")
        )
      )
  }
}

if (!"lux_plot_f" %in% names(pred_exp1)) {
  if ("lux_exposure" %in% names(pred_exp1)) {
    pred_exp1 <- pred_exp1 |>
      dplyr::mutate(
        lux_plot_f = factor(
          paste0(lux_exposure, " lux"),
          levels = c("0 lux", "4 lux", "70 lux", "117 lux")
        )
      )
  }
}

check_required_columns(
  raw_exp1,
  c("motility_ratio", "lux_plot_f"),
  "Experiment 1 raw data"
)

check_required_columns(
  pred_exp1,
  c("lux_plot_f", "fit_response", "conf_low_response", "conf_high_response"),
  "Experiment 1 predictions"
)

warn_if_zoom_clips(
  values_low = c(raw_exp1$motility_ratio, pred_exp1$conf_low_response),
  values_high = c(raw_exp1$motility_ratio, pred_exp1$conf_high_response),
  zoom_limits = zoom_y_limits,
  label = "Figure 4"
)

p_fig4_zoom <- ggplot2::ggplot(
  raw_exp1,
  ggplot2::aes(x = lux_plot_f, y = motility_ratio)
) +
  ggplot2::geom_point(
    colour = kelp_raw_point_colour,
    alpha = 0.45,
    size = 2.2,
    position = ggplot2::position_jitter(width = 0.08, height = 0)
  ) +
  ggplot2::geom_errorbar(
    data = pred_exp1,
    ggplot2::aes(
      x = lux_plot_f,
      y = fit_response,
      ymin = conf_low_response,
      ymax = conf_high_response,
      colour = lux_plot_f
    ),
    width = 0.08,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = pred_exp1,
    ggplot2::aes(
      x = lux_plot_f,
      y = fit_response,
      colour = lux_plot_f
    ),
    size = 3.2,
    inherit.aes = FALSE
  ) +
  scale_colour_kelp_exp1_light(drop = FALSE) +
  ggplot2::coord_cartesian(ylim = zoom_y_limits) +
  ggplot2::labs(
    title = "Experiment 1: Light-only gradient",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Irradiance treatment",
    y = "Motile fraction",
    colour = "Irradiance treatment",
    caption = paste(
      "Points show raw well-level observations.",
      "Coloured points and error bars show model-predicted motile fraction with 95% confidence intervals.",
      "Display window zoomed for comparison; statistical outputs unchanged."
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5)
  )

fig4_manifest <- save_plot_multi_local(
  plot = p_fig4_zoom,
  file_stem = "Fig4_exp1_light_day4_model_predictions_zoomed",
  width = fig_width_std,
  height = fig_height_std,
  dpi = fig_dpi
)

cat("Figure 4 zoomed saved.\n")
print(fig4_manifest)
cat("\n")

# ---------------------------------------------------------
# 10. Figure 5 - Experiment 2 zoomed
# ---------------------------------------------------------

raw_exp2 <- day4_data |>
  dplyr::filter(experiment_num == 10.2)

# Flexible label prep
if (!"ntu_plot_f" %in% names(raw_exp2)) {
  if ("ntu" %in% names(raw_exp2)) {
    raw_exp2 <- raw_exp2 |>
      dplyr::mutate(
        ntu_plot_f = factor(
          paste0(ntu, " NTU"),
          levels = c("25 NTU", "100 NTU", "400 NTU")
        )
      )
  }
}

if (!"particle_type_plot" %in% names(raw_exp2)) {
  if ("particle_type" %in% names(raw_exp2)) {
    raw_exp2 <- raw_exp2 |>
      dplyr::mutate(
        particle_type_plot = dplyr::case_when(
          stringr::str_to_upper(as.character(particle_type)) == "SAND" ~ "Sand",
          stringr::str_to_upper(as.character(particle_type)) == "KAOLINITE" ~ "Kaolinite",
          stringr::str_to_upper(as.character(particle_type)) == "PEAT" ~ "Peat",
          TRUE ~ as.character(particle_type)
        )
      )
  }
}

if (!"ntu_plot_f" %in% names(pred_exp2)) {
  if ("ntu" %in% names(pred_exp2)) {
    pred_exp2 <- pred_exp2 |>
      dplyr::mutate(
        ntu_plot_f = factor(
          paste0(ntu, " NTU"),
          levels = c("25 NTU", "100 NTU", "400 NTU")
        )
      )
  }
}

check_required_columns(
  raw_exp2,
  c("motility_ratio", "ntu_plot_f", "particle_type_plot"),
  "Experiment 2 raw data"
)

check_required_columns(
  pred_exp2,
  c("ntu_plot_f", "particle_type_plot", "fit_response", "conf_low_response", "conf_high_response"),
  "Experiment 2 predictions"
)

warn_if_zoom_clips(
  values_low = c(raw_exp2$motility_ratio, pred_exp2$conf_low_response),
  values_high = c(raw_exp2$motility_ratio, pred_exp2$conf_high_response),
  zoom_limits = zoom_y_limits,
  label = "Figure 5"
)

pd_exp2 <- ggplot2::position_dodge(width = 0.38)
pjd_exp2 <- ggplot2::position_jitterdodge(
  jitter.width = 0.08,
  dodge.width = 0.38,
  jitter.height = 0
)

p_fig5_zoom <- ggplot2::ggplot(
  raw_exp2,
  ggplot2::aes(
    x = ntu_plot_f,
    y = motility_ratio,
    colour = particle_type_plot
  )
) +
  ggplot2::geom_errorbar(
    data = pred_exp2,
    ggplot2::aes(
      x = ntu_plot_f,
      ymin = conf_low_response,
      ymax = conf_high_response,
      colour = particle_type_plot,
      group = particle_type_plot
    ),
    width = 0.08,
    linewidth = 0.8,
    position = pd_exp2,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = pred_exp2,
    ggplot2::aes(
      x = ntu_plot_f,
      y = fit_response,
      colour = particle_type_plot,
      group = particle_type_plot
    ),
    size = 3.0,
    position = pd_exp2,
    inherit.aes = FALSE
  ) +
  scale_colour_kelp_exp2_particle(drop = FALSE) +
  ggplot2::coord_cartesian(ylim = zoom_y_limits) +
  ggplot2::labs(
    title = "Experiment 2: Defined particle concentration series",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Turbidity treatment",
    y = "Motile fraction",
    colour = "Particle type",
    caption = paste(
      "Points show raw well-level observations.",
      "Model predictions are shown for each defined particle type at 25, 100, and 400 NTU.",
      "Display window zoomed for comparison; statistical outputs unchanged."
    )
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5)
  )

fig5_manifest <- save_plot_multi_local(
  plot = p_fig5_zoom,
  file_stem = "Fig5_exp2_defined_particles_day4_model_predictions_zoomed",
  width = fig_width_std,
  height = fig_height_std,
  dpi = fig_dpi
)

cat("Figure 5 zoomed saved.\n")
print(fig5_manifest)
cat("\n")

# ---------------------------------------------------------
# 11. Figure 6 - Experiment 3 zoomed
# ---------------------------------------------------------

raw_exp3 <- day4_data |>
  dplyr::filter(experiment_num == 11.2)

# Prepare brake-wear load labels for raw data
if (!"bw_load_f" %in% names(raw_exp3)) {
  if ("bw_ug_l" %in% names(raw_exp3)) {
    raw_exp3 <- raw_exp3 |>
      dplyr::mutate(
        bw_load_f = factor(
          paste0(bw_ug_l, " ug/L"),
          levels = c("4.5 ug/L", "45 ug/L", "450 ug/L")
        )
      )
  }
}

# Prepare brake-wear load labels for prediction data
if (!"bw_load_f" %in% names(pred_exp3)) {
  if ("bw_ug_l" %in% names(pred_exp3)) {
    pred_exp3 <- pred_exp3 |>
      dplyr::mutate(
        bw_load_f = factor(
          paste0(bw_ug_l, " ug/L"),
          levels = c("4.5 ug/L", "45 ug/L", "450 ug/L")
        )
      )
  }
}

# Prepare size/fraction labels for raw data
if (!"particle_size_plot" %in% names(raw_exp3)) {
  
  size_source_col <- dplyr::case_when(
    "particle_size_class" %in% names(raw_exp3) ~ "particle_size_class",
    "size_class_plot" %in% names(raw_exp3) ~ "size_class_plot",
    "size_class" %in% names(raw_exp3) ~ "size_class",
    "brake_size_treatment" %in% names(raw_exp3) ~ "brake_size_treatment",
    "toxin_exposure" %in% names(raw_exp3) ~ "toxin_exposure",
    TRUE ~ NA_character_
  )
  
  if (is.na(size_source_col)) {
    stop(
      paste0(
        "Experiment 3 raw data is missing a usable brake-wear size/fraction column. ",
        "Expected one of particle_size_class, size_class_plot, size_class, ",
        "brake_size_treatment, or toxin_exposure."
      ),
      call. = FALSE
    )
  }
  
  raw_exp3 <- raw_exp3 |>
    dplyr::mutate(
      size_source = .data[[size_source_col]],
      size_source_clean = stringr::str_to_upper(
        stringr::str_squish(as.character(size_source))
      ),
      particle_size_plot = dplyr::case_when(
        size_source_clean %in% c(
          "BWC",
          "COARSE",
          "BRAKE WEAR COARSE",
          "BRAKE_WEAR_COARSE",
          "BRAKE COARSE"
        ) ~ "Coarse",
        size_source_clean %in% c(
          "BWF",
          "FINE",
          "BRAKE WEAR FINE",
          "BRAKE_WEAR_FINE",
          "BRAKE FINE"
        ) ~ "Fine",
        stringr::str_detect(size_source_clean, "BWC|COARSE") ~ "Coarse",
        stringr::str_detect(size_source_clean, "BWF|FINE") ~ "Fine",
        TRUE ~ NA_character_
      ),
      particle_size_plot = factor(
        particle_size_plot,
        levels = c("Coarse", "Fine")
      )
    )
  
  cat("Experiment 3 raw size/fraction source column:", size_source_col, "\n")
}

# Prepare size/fraction labels for prediction data
if (!"particle_size_plot" %in% names(pred_exp3)) {
  
  pred_size_source_col <- dplyr::case_when(
    "particle_size_class" %in% names(pred_exp3) ~ "particle_size_class",
    "size_class_plot" %in% names(pred_exp3) ~ "size_class_plot",
    "size_class" %in% names(pred_exp3) ~ "size_class",
    "brake_size_treatment" %in% names(pred_exp3) ~ "brake_size_treatment",
    TRUE ~ NA_character_
  )
  
  if (is.na(pred_size_source_col)) {
    stop(
      paste0(
        "Experiment 3 prediction data is missing a usable brake-wear ",
        "size/fraction column. Expected one of particle_size_class, ",
        "size_class_plot, size_class, or brake_size_treatment."
      ),
      call. = FALSE
    )
  }
  
  pred_exp3 <- pred_exp3 |>
    dplyr::mutate(
      pred_size_source = .data[[pred_size_source_col]],
      pred_size_source_clean = stringr::str_to_upper(
        stringr::str_squish(as.character(pred_size_source))
      ),
      particle_size_plot = dplyr::case_when(
        pred_size_source_clean %in% c("BWC", "COARSE") ~ "Coarse",
        pred_size_source_clean %in% c("BWF", "FINE") ~ "Fine",
        stringr::str_detect(pred_size_source_clean, "BWC|COARSE") ~ "Coarse",
        stringr::str_detect(pred_size_source_clean, "BWF|FINE") ~ "Fine",
        TRUE ~ as.character(pred_size_source)
      )
    )
  
  cat("Experiment 3 prediction size/fraction source column:", pred_size_source_col, "\n")
}

# Standardise prediction size/fraction labels and factor order
pred_exp3 <- pred_exp3 |>
  dplyr::mutate(
    particle_size_plot = as.character(particle_size_plot),
    particle_size_plot = dplyr::case_when(
      particle_size_plot %in% c("Coarse", "coarse", "BWC") ~ "Coarse",
      particle_size_plot %in% c("Fine", "fine", "BWF") ~ "Fine",
      TRUE ~ particle_size_plot
    ),
    particle_size_plot = factor(
      particle_size_plot,
      levels = c("Coarse", "Fine")
    )
  )

# Check required columns before plotting
check_required_columns(
  raw_exp3,
  c("motility_ratio", "bw_load_f", "particle_size_plot"),
  "Experiment 3 raw data"
)

check_required_columns(
  pred_exp3,
  c(
    "bw_load_f",
    "particle_size_plot",
    "fit_response",
    "conf_low_response",
    "conf_high_response"
  ),
  "Experiment 3 predictions"
)

# Warn if the zoom window excludes visible values
warn_if_zoom_clips(
  values_low = c(raw_exp3$motility_ratio, pred_exp3$conf_low_response),
  values_high = c(raw_exp3$motility_ratio, pred_exp3$conf_high_response),
  zoom_limits = zoom_y_limits,
  label = "Figure 6"
)

pd_exp3 <- ggplot2::position_dodge(width = 0.35)

pjd_exp3 <- ggplot2::position_jitterdodge(
  jitter.width = 0.08,
  dodge.width = 0.35,
  jitter.height = 0
)

p_fig6_zoom <- ggplot2::ggplot(
  raw_exp3,
  ggplot2::aes(
    x = bw_load_f,
    y = motility_ratio,
    colour = particle_size_plot
  )
) +
  ggplot2::geom_point(
    alpha = 0.35,
    size = 2.2,
    position = pjd_exp3
  ) +
  ggplot2::geom_errorbar(
    data = pred_exp3,
    ggplot2::aes(
      x = bw_load_f,
      y = fit_response,
      ymin = conf_low_response,
      ymax = conf_high_response,
      colour = particle_size_plot,
      group = particle_size_plot
    ),
    width = 0.08,
    linewidth = 0.8,
    position = pd_exp3,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = pred_exp3,
    ggplot2::aes(
      x = bw_load_f,
      y = fit_response,
      colour = particle_size_plot,
      group = particle_size_plot
    ),
    size = 3.0,
    position = pd_exp3,
    inherit.aes = FALSE
  ) +
  scale_colour_kelp_exp3_brake(drop = FALSE) +
  ggplot2::coord_cartesian(ylim = zoom_y_limits) +
  ggplot2::labs(
    title = "Experiment 3: Brake-wear size comparison",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "Brake-wear load",
    y = "Motile fraction",
    colour = "Brake-wear fraction",
    caption = paste(
      "Points show raw well-level observations.",
      "Model predictions are shown for coarse and fine brake-wear fractions across mass loads.",
      "Display window zoomed for comparison; statistical outputs unchanged."
    )
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5)
  )

fig6_manifest <- save_plot_multi_local(
  plot = p_fig6_zoom,
  file_stem = "Fig6_exp3_brake_size_day4_model_predictions_zoomed",
  width = fig_width_std,
  height = fig_height_std,
  dpi = fig_dpi
)

cat("Figure 6 zoomed saved.\n")
print(fig6_manifest)
cat("\n")
# ---------------------------------------------------------
# 12. Figure 7 - Experiment 4 zoomed
# ---------------------------------------------------------

raw_exp4 <- day4_data |>
  dplyr::filter(experiment_num == 8.2)

if (!"log10_ntu_plus1" %in% names(raw_exp4)) {
  if ("ntu" %in% names(raw_exp4)) {
    raw_exp4 <- raw_exp4 |>
      dplyr::mutate(log10_ntu_plus1 = log10(ntu + 1))
  }
}

if (!"log10_ntu_plus1" %in% names(pred_exp4)) {
  if ("ntu" %in% names(pred_exp4)) {
    pred_exp4 <- pred_exp4 |>
      dplyr::mutate(log10_ntu_plus1 = log10(ntu + 1))
  }
}

check_required_columns(
  raw_exp4,
  c("motility_ratio", "log10_ntu_plus1"),
  "Experiment 4 raw data"
)

check_required_columns(
  pred_exp4,
  c("log10_ntu_plus1", "fit_response", "conf_low_response", "conf_high_response"),
  "Experiment 4 predictions"
)

warn_if_zoom_clips(
  values_low = c(raw_exp4$motility_ratio, pred_exp4$conf_low_response),
  values_high = c(raw_exp4$motility_ratio, pred_exp4$conf_high_response),
  zoom_limits = zoom_y_limits,
  label = "Figure 7"
)

spm_colour <- unname(kelp_exp4_spm_family[[1]])

p_fig7_zoom <- ggplot2::ggplot(
  raw_exp4,
  ggplot2::aes(x = log10_ntu_plus1, y = motility_ratio)
) +
  ggplot2::geom_point(
    colour = kelp_raw_point_colour,
    alpha = 0.45,
    size = 2.3
  ) +
  ggplot2::geom_ribbon(
    data = pred_exp4,
    ggplot2::aes(
      x = log10_ntu_plus1,
      ymin = conf_low_response,
      ymax = conf_high_response
    ),
    fill = spm_colour,
    alpha = 0.18,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_line(
    data = pred_exp4,
    ggplot2::aes(
      x = log10_ntu_plus1,
      y = fit_response
    ),
    colour = spm_colour,
    linewidth = 1.0,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_cartesian(ylim = zoom_y_limits) +
  ggplot2::labs(
    title = "Experiment 4: Field-derived SPM gradient",
    subtitle = "Day 4 model-predicted motile fraction with 95% confidence intervals",
    x = "log10(NTU + 1)",
    y = "Motile fraction",
    caption = paste(
      "Points show raw well-level observations.",
      "Line and ribbon show saved model predictions with 95% confidence intervals.",
      "Display window zoomed for comparison; statistical outputs unchanged."
    )
  ) +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5)
  )

fig7_manifest <- save_plot_multi_local(
  plot = p_fig7_zoom,
  file_stem = "Fig7_exp4_field_spm_day4_model_predictions_zoomed",
  width = fig_width_std,
  height = fig_height_std,
  dpi = fig_dpi
)

cat("Figure 7 zoomed saved.\n")
print(fig7_manifest)
cat("\n")

# ---------------------------------------------------------
# 13. Figure 8 - zoomed synthesis forest plot
# ---------------------------------------------------------

check_required_columns(
  fig8_combined,
  c(
    "experiment_block",
    "display_label",
    "percent_change",
    "conf_low_percent",
    "conf_high_percent"
  ),
  "Figure 8 combined effects"
)

# Rebuild ordering safely
fig8_data <- fig8_combined |>
  dplyr::mutate(
    experiment_block = factor(
      experiment_block,
      levels = c(
        "Experiment 1: light-only",
        "Experiment 2: defined particles",
        "Experiment 3: brake-wear size",
        "Experiment 4: field-derived SPM"
      )
    ),
    display_label = as.character(display_label),
    ci_crosses_zero = if ("ci_crosses_zero" %in% names(fig8_combined)) {
      as.logical(ci_crosses_zero)
    } else {
      conf_low_percent <= 0 & conf_high_percent >= 0
    },
    within_experiment_order = dplyr::case_when(
      experiment_block == "Experiment 1: light-only" & display_label == "0 lux" ~ 1,
      experiment_block == "Experiment 1: light-only" & display_label == "4 lux" ~ 2,
      experiment_block == "Experiment 1: light-only" & display_label == "70 lux" ~ 3,
      
      experiment_block == "Experiment 2: defined particles" & display_label == "Kaolinite" ~ 1,
      experiment_block == "Experiment 2: defined particles" & display_label == "Peat" ~ 2,
      
      experiment_block == "Experiment 3: brake-wear size" & display_label == "Coarse 45 ug/L" ~ 1,
      experiment_block == "Experiment 3: brake-wear size" & display_label == "Coarse 450 ug/L" ~ 2,
      experiment_block == "Experiment 3: brake-wear size" & display_label == "Fine 4.5 ug/L" ~ 3,
      experiment_block == "Experiment 3: brake-wear size" & display_label == "Fine 45 ug/L" ~ 4,
      experiment_block == "Experiment 3: brake-wear size" & display_label == "Fine 450 ug/L" ~ 5,
      
      experiment_block == "Experiment 4: field-derived SPM" & display_label == "25 NTU" ~ 1,
      experiment_block == "Experiment 4: field-derived SPM" & display_label == "100 NTU" ~ 2,
      experiment_block == "Experiment 4: field-derived SPM" & display_label == "400 NTU" ~ 3,
      experiment_block == "Experiment 4: field-derived SPM" & display_label == "500 NTU" ~ 4,
      experiment_block == "Experiment 4: field-derived SPM" & display_label == "4144 NTU (extreme)" ~ 5,
      
      TRUE ~ 999
    )
  ) |>
  dplyr::arrange(experiment_block, within_experiment_order) |>
  dplyr::group_by(experiment_block) |>
  dplyr::mutate(display_order = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    display_label_unique = paste0(
      display_label,
      "___",
      as.integer(experiment_block),
      "_",
      sprintf("%02d", display_order)
    ),
    display_label_unique = factor(
      display_label_unique,
      levels = rev(unique(display_label_unique))
    )
  )

warn_if_zoom_clips(
  values_low = fig8_data$conf_low_percent,
  values_high = fig8_data$conf_high_percent,
  zoom_limits = zoom_x_limits_fig8,
  label = "Figure 8"
)

clean_y_labels <- function(x) {
  stringr::str_replace(as.character(x), "___.*$", "")
}

p_fig8_zoom <- ggplot2::ggplot(
  fig8_data,
  ggplot2::aes(
    x = percent_change,
    y = display_label_unique
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    colour = kelp_reference_line_colour
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = conf_low_percent,
      xmax = conf_high_percent,
      colour = experiment_block
    ),
    height = 0.18,
    linewidth = 0.75,
    alpha = 0.95
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = experiment_block,
      shape = ci_crosses_zero
    ),
    size = 2.9
  ) +
  ggplot2::facet_grid(
    experiment_block ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_y_discrete(labels = clean_y_labels) +
  scale_colour_kelp_experiment(drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 16),
    labels = c(`TRUE` = "CI crosses 0", `FALSE` = "CI excludes 0"),
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(xlim = zoom_x_limits_fig8) +
  ggplot2::labs(
    title = "Figure 8: Relative Day 4 effects across experiments",
    subtitle = "Predicted proportional change in motile fraction relative to experiment-specific reference conditions",
    x = "Predicted change in motile fraction relative to reference (%)",
    y = NULL,
    colour = "Experiment",
    shape = "95% CI",
    caption = paste(
      "Negative values indicate reduced predicted motile fraction; positive values indicate increased predicted motile fraction.",
      "Reference conditions: 117 lux; sand at 25 NTU; BWC at 4.5 ug/L; and lowest observed field-SPM NTU.",
      "Display window zoomed for comparison; statistical outputs unchanged."
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      hjust = 0
    ),
    panel.spacing.y = grid::unit(0.7, "lines"),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8.5)
  )

fig8_manifest <- save_plot_multi_local(
  plot = p_fig8_zoom,
  file_stem = "Fig8_day4_relative_effects_forest_plot_zoomed",
  width = fig_width_f8,
  height = fig_height_f8,
  dpi = fig_dpi
)

cat("Figure 8 zoomed saved.\n")
print(fig8_manifest)
cat("\n")

# ---------------------------------------------------------
# 14. Manifest
# ---------------------------------------------------------

zoom_manifest <- dplyr::bind_rows(
  fig4_manifest |>
    dplyr::mutate(
      figure_role = "Day 4 main figure",
      zoom_axis = "y",
      zoom_limits = paste(zoom_y_limits, collapse = " to "),
      source_file = file_pred_exp1
    ),
  fig5_manifest |>
    dplyr::mutate(
      figure_role = "Day 4 main figure",
      zoom_axis = "y",
      zoom_limits = paste(zoom_y_limits, collapse = " to "),
      source_file = file_pred_exp2
    ),
  fig6_manifest |>
    dplyr::mutate(
      figure_role = "Day 4 main figure",
      zoom_axis = "y",
      zoom_limits = paste(zoom_y_limits, collapse = " to "),
      source_file = file_pred_exp3
    ),
  fig7_manifest |>
    dplyr::mutate(
      figure_role = "Day 4 main figure",
      zoom_axis = "y",
      zoom_limits = paste(zoom_y_limits, collapse = " to "),
      source_file = file_pred_exp4
    ),
  fig8_manifest |>
    dplyr::mutate(
      figure_role = "Day 4 synthesis figure",
      zoom_axis = "x",
      zoom_limits = paste(zoom_x_limits_fig8, collapse = " to "),
      source_file = file_fig8_combined
    )
)

readr::write_csv(zoom_manifest, file_manifest)

cat("Zoomed figure manifest written to:\n")
cat(file_manifest, "\n\n")

print(zoom_manifest)

# ---------------------------------------------------------
# 15. Log
# ---------------------------------------------------------

sink(file_log)

cat("Zoomed Day 4 figures log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input files:\n")
cat("Day 4 dataset:", if (file.exists(file_day4_rds)) file_day4_rds else file_day4_csv, "\n")
cat("Exp1 predictions:", file_pred_exp1, "\n")
cat("Exp2 predictions:", file_pred_exp2, "\n")
cat("Exp3 predictions:", file_pred_exp3, "\n")
cat("Exp4 predictions:", file_pred_exp4, "\n")
cat("Figure 8 combined effects:", file_fig8_combined, "\n\n")

cat("Zoom settings:\n")
cat("Figures 4-7 y-axis:", paste(zoom_y_limits, collapse = " to "), "\n")
cat("Figure 8 x-axis:", paste(zoom_x_limits_fig8, collapse = " to "), "\n\n")

cat("Rows loaded:\n")
cat("Day 4 dataset:", nrow(day4_data), "\n")
cat("Exp1 predictions:", nrow(pred_exp1), "\n")
cat("Exp2 predictions:", nrow(pred_exp2), "\n")
cat("Exp3 predictions:", nrow(pred_exp3), "\n")
cat("Exp4 predictions:", nrow(pred_exp4), "\n")
cat("Figure 8 combined:", nrow(fig8_data), "\n\n")

cat("Manifest:\n")
print(zoom_manifest)

sink()

cat("Log written to:\n")
cat(file_log, "\n\n")

cat("SCRIPT 14 COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
