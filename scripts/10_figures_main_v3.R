# =========================================================
# Script title: 10_figures_main_v3.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Generate the main publication figures for
#          Manuscript 4 from the derived dataset and saved
#          experiment-specific model outputs.
# Inputs:
# - data_derived/ms4_analysis_derived.csv or .rds
# - outputs/models/exp1_light/exp1_light_preferred_model_predictions.csv
# - outputs/models/exp2_defined_particles/exp2_defined_particles_preferred_model_predictions.csv
# - outputs/models/exp3_brake_size/exp3_brake_size_preferred_model_predictions.csv
# - outputs/models/exp4_field_spm/exp4_field_spm_preferred_model_predictions.csv
# - outputs/models/exp4_field_spm/exp4_field_spm_day4_preferred_model_predictions.csv
# Outputs:
# - outputs/figures/models_exp1/Fig2_light_model_predictions.{pdf,png,tiff}
# - outputs/figures/models_exp2/Fig3_defined_particles_model_predictions.{pdf,png,tiff}
# - outputs/figures/models_exp3/Fig4_brake_size_model_predictions.{pdf,png,tiff}
# - outputs/figures/models_exp4/Fig5_field_spm_model_predictions.{pdf,png,tiff}
# - outputs/figures/models_exp4/FigSx_field_spm_day4_sensitivity.{pdf,png,tiff}
# - paired caption markdown files
# - outputs/logs/10_figures_main_log_*.txt
# Date created: 26 March 2026
# Last updated: 31 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 06-09 model scripts before this script.
# - Manuscript figures follow the publication plan:
#   Fig 2 Light-only gradient
#   Fig 3 Defined particle concentration series
#   Fig 4 Brake-wear size comparison
#   Fig 5 Field-derived SPM gradient
# - Experiment 4 Day 4 sensitivity is exported separately
#   for optional manual assembly into the final figure.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 10: FIGURES MAIN\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check setup objects
# ---------------------------------------------------------

required_objects <- c(
  "project_root",
  "dir_data_derived",
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
  "helpers_labels.R"
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
# 3. Package checks
# ---------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "ggplot2", "stringr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "The following package(s) are required but not installed:\n- ",
      paste(missing_pkgs, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------
# 4. Small local helpers
# ---------------------------------------------------------

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
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

write_caption_md_local <- function(figure_name, caption_text, subdir) {
  fig_dir <- file.path(project_root, "outputs", "figures", subdir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
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

save_plot_multi <- function(plot,
                            figure_name,
                            subdir,
                            width = 8,
                            height = 6,
                            dpi = 600,
                            caption_text = NULL) {
  fig_dir <- file.path(project_root, "outputs", "figures", subdir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  file_png <- file.path(fig_dir, paste0(figure_name, ".png"))
  file_tiff <- file.path(fig_dir, paste0(figure_name, ".tiff"))
  file_pdf <- file.path(fig_dir, paste0(figure_name, ".pdf"))
  
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
  
  ggplot2::ggsave(
    filename = file_pdf,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  caption_file <- NULL
  if (!is.null(caption_text) && nzchar(caption_text)) {
    caption_file <- if (exists("write_figure_caption_md", mode = "function", inherits = TRUE)) {
      write_figure_caption_md(
        figure_name = figure_name,
        caption_text = caption_text,
        subdir = subdir
      )
    } else {
      write_caption_md_local(
        figure_name = figure_name,
        caption_text = caption_text,
        subdir = subdir
      )
    }
  }
  
  tibble::tibble(
    figure_name = figure_name,
    subdir = subdir,
    file_png = file_png,
    file_tiff = file_tiff,
    file_pdf = file_pdf,
    caption_file = caption_file
  )
}

standardise_upper_trim <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
}

# ---------------------------------------------------------
# 5. Load derived dataset
# ---------------------------------------------------------

file_input_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")

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
      file_input_csv
    ),
    call. = FALSE
  )
}

cat("Derived dataset loaded from:\n")
cat(input_source_used, "\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

# ---------------------------------------------------------
# 6. Restrict to Manuscript 4 block and derive plotting vars
# ---------------------------------------------------------

ms4_experiment_nums <- c(8.2, 9.2, 10.2, 11.2)

dat <- dat |>
  dplyr::mutate(
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    ntu = if ("ntu" %in% names(.)) suppressWarnings(as.numeric(ntu)) else NA_real_,
    bw_ug_l = if ("bw_ug_l" %in% names(.)) suppressWarnings(as.numeric(bw_ug_l)) else NA_real_,
    lux_exposure = if ("lux_exposure" %in% names(.)) suppressWarnings(as.numeric(lux_exposure)) else NA_real_,
    days_from_start = if ("days_from_start" %in% names(.)) suppressWarnings(as.integer(days_from_start)) else NA_integer_
  ) |>
  dplyr::filter(
    !is.na(experiment_num),
    experiment_num %in% ms4_experiment_nums
  ) |>
  dplyr::mutate(
    total_cells = mobile_cell_count + stationary_cell_count,
    motility_ratio = dplyr::if_else(
      total_cells > 0,
      mobile_cell_count / total_cells,
      NA_real_
    ),
    log10_ntu_plus1 = ifelse(!is.na(ntu), log10(ntu + 1), NA_real_),
    log10_bw_ug_l_plus1 = ifelse(!is.na(bw_ug_l), log10(bw_ug_l + 1), NA_real_),
    day_label = ifelse(!is.na(days_from_start), paste0("Day ", days_from_start), "Single day")
  )

cat("Restricted to Manuscript 4 block.\n\n")

# ---------------------------------------------------------
# 7. Load saved model predictions
# ---------------------------------------------------------

pred_exp1 <- read_if_exists(file.path(
  project_root, "outputs", "models", "exp1_light",
  "exp1_light_preferred_model_predictions.csv"
))

pred_exp2 <- read_if_exists(file.path(
  project_root, "outputs", "models", "exp2_defined_particles",
  "exp2_defined_particles_preferred_model_predictions.csv"
))

pred_exp3 <- read_if_exists(file.path(
  project_root, "outputs", "models", "exp3_brake_size",
  "exp3_brake_size_preferred_model_predictions.csv"
))

pred_exp4 <- read_if_exists(file.path(
  project_root, "outputs", "models", "exp4_field_spm",
  "exp4_field_spm_preferred_model_predictions.csv"
))

pred_exp4_day4 <- read_if_exists(file.path(
  project_root, "outputs", "models", "exp4_field_spm",
  "exp4_field_spm_day4_preferred_model_predictions.csv"
))

# ---------------------------------------------------------
# 8. Define palettes
# ---------------------------------------------------------

particle_palette_local <- c(
  "CONTROL" = "#4D4D4D",
  "SAND" = "#C2B280",
  "KAOLINITE" = "#A6CEE3",
  "PEAT" = "#8B4513",
  "BWC" = "#D55E00",
  "BWF" = "#CC79A7",
  "SPM" = "#0072B2"
)

light_palette_local <- c(
  "0" = "#f7f7f7",
  "4" = "#cccccc",
  "70" = "#969696",
  "117" = "#525252"
)

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

manifest_figures <- list()

# ---------------------------------------------------------
# 9. Figure 2: Experiment 1
# ---------------------------------------------------------

exp1_raw <- dat |>
  dplyr::filter(experiment_num == 9.2)

if (!is.null(pred_exp1) && nrow(exp1_raw) > 0) {
  pred_exp1 <- pred_exp1 |>
    dplyr::mutate(
      days_from_start = suppressWarnings(as.integer(days_from_start)),
      day_label = paste0("Day ", days_from_start)
    )
  
  exp1_raw <- exp1_raw |>
    dplyr::mutate(
      lux_exposure_f = factor(as.character(lux_exposure), levels = c("0", "4", "70", "117")),
      day_label = paste0("Day ", days_from_start)
    )
  
  p_fig2 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = exp1_raw,
      ggplot2::aes(
        x = lux_exposure,
        y = motility_ratio,
        colour = lux_exposure_f
      ),
      alpha = 0.35,
      size = 1.7
    ) +
    ggplot2::geom_line(
      data = pred_exp1,
      ggplot2::aes(
        x = lux_exposure,
        y = fit_response
      ),
      linewidth = 0.9,
      colour = "#404040"
    ) +
    ggplot2::facet_wrap(~ day_label) +
    ggplot2::scale_colour_manual(values = light_palette_local, drop = FALSE) +
    ggplot2::labs(
      title = "Light-only irradiance gradient",
      subtitle = "Model-predicted motile fraction versus irradiance",
      x = "Irradiance (lux)",
      y = "Motile fraction",
      colour = "Lux",
      caption = "Points show raw observations and lines show model predictions."
    )
  
  manifest_figures[[length(manifest_figures) + 1]] <- save_plot_multi(
    plot = p_fig2,
    figure_name = "Fig2_light_model_predictions",
    subdir = "models_exp1",
    width = 8,
    height = 6.5,
    caption_text = paste(
      "Figure 2. Light-only irradiance gradient.",
      "Model-predicted motile fraction as a function of irradiance, facetted by day.",
      "Points show raw observations and lines show model predictions."
    )
  )
  
  cat("Figure 2 saved.\n")
}

# ---------------------------------------------------------
# 10. Figure 3: Experiment 2
# ---------------------------------------------------------

exp2_raw <- dat |>
  dplyr::filter(experiment_num == 10.2)

if (!is.null(pred_exp2) && nrow(exp2_raw) > 0) {
  exp2_raw <- exp2_raw |>
    dplyr::mutate(
      particle_type_std = standardise_upper_trim(particle_type),
      particle_type_std = dplyr::case_when(
        particle_type_std %in% c("CONTROL") ~ "CONTROL",
        particle_type_std %in% c("SAND") ~ "SAND",
        particle_type_std %in% c("KAOLINITE", "KAOLIN") ~ "KAOLINITE",
        particle_type_std %in% c("PEAT") ~ "PEAT",
        TRUE ~ particle_type_std
      ),
      day_label = paste0("Day ", days_from_start)
    )
  
  pred_exp2 <- pred_exp2 |>
    dplyr::mutate(
      particle_type = standardise_upper_trim(particle_type),
      particle_type = dplyr::case_when(
        particle_type %in% c("CONTROL") ~ "CONTROL",
        particle_type %in% c("SAND") ~ "SAND",
        particle_type %in% c("KAOLINITE", "KAOLIN") ~ "KAOLINITE",
        particle_type %in% c("PEAT") ~ "PEAT",
        TRUE ~ particle_type
      ),
      days_from_start = suppressWarnings(as.integer(days_from_start)),
      day_label = paste0("Day ", days_from_start)
    )
  
  p_fig3 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = exp2_raw,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = motility_ratio,
        colour = particle_type_std
      ),
      alpha = 0.28,
      size = 1.6
    ) +
    ggplot2::geom_line(
      data = pred_exp2,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = fit_response,
        colour = particle_type
      ),
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(~ day_label) +
    ggplot2::scale_colour_manual(values = particle_palette_local, drop = FALSE) +
    ggplot2::labs(
      title = "Defined particle concentration series",
      subtitle = "Model-predicted motile fraction versus log10(NTU + 1), by particle type",
      x = "log10(NTU + 1)",
      y = "Motile fraction",
      colour = "Particle type",
      caption = "Points show raw observations and lines show model predictions."
    )
  
  manifest_figures[[length(manifest_figures) + 1]] <- save_plot_multi(
    plot = p_fig3,
    figure_name = "Fig3_defined_particles_model_predictions",
    subdir = "models_exp2",
    width = 8.5,
    height = 6.5,
    caption_text = paste(
      "Figure 3. Defined particle concentration series.",
      "Model-predicted motile fraction as a function of log10(NTU + 1), coloured by particle type and facetted by day.",
      "Points show raw observations and lines show model predictions."
    )
  )
  
  cat("Figure 3 saved.\n")
}

# ---------------------------------------------------------
# 11. Figure 4: Experiment 3
# ---------------------------------------------------------

exp3_raw <- dat |>
  dplyr::filter(experiment_num == 11.2) |>
  dplyr::mutate(
    particle_size_class = standardise_upper_trim(
      dplyr::coalesce(
        if ("size_class_plot" %in% names(.)) size_class_plot else NA_character_,
        if ("size_class" %in% names(.)) size_class else NA_character_
      )
    ),
    particle_size_class = dplyr::case_when(
      particle_size_class %in% c("COARSE", "BWC", "BRAKE WEAR COARSE", "BRAKE_WEAR_COARSE") ~ "BWC",
      particle_size_class %in% c("FINE", "BWF", "BRAKE WEAR FINE", "BRAKE_WEAR_FINE") ~ "BWF",
      TRUE ~ NA_character_
    )
  )

if (!is.null(pred_exp3) && nrow(exp3_raw) > 0) {
  pred_exp3 <- pred_exp3 |>
    dplyr::mutate(
      particle_size_class = standardise_upper_trim(particle_size_class)
    )
  
  p_fig4 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = exp3_raw,
      ggplot2::aes(
        x = log10_bw_ug_l_plus1,
        y = motility_ratio,
        colour = particle_size_class
      ),
      alpha = 0.30,
      size = 1.8
    ) +
    ggplot2::geom_line(
      data = pred_exp3,
      ggplot2::aes(
        x = log10_bw_ug_l_plus1,
        y = fit_response,
        colour = particle_size_class
      ),
      linewidth = 0.9
    ) +
    ggplot2::scale_colour_manual(values = particle_palette_local, drop = FALSE) +
    ggplot2::labs(
      title = "Brake-wear size comparison under equal mass loading",
      subtitle = "Model-predicted motile fraction versus log10(brake-wear concentration + 1)",
      x = "log10(brake-wear concentration [ug/L] + 1)",
      y = "Motile fraction",
      colour = "Particle class",
      caption = "Points show raw observations and lines show model predictions."
    )
  
  manifest_figures[[length(manifest_figures) + 1]] <- save_plot_multi(
    plot = p_fig4,
    figure_name = "Fig4_brake_size_model_predictions",
    subdir = "models_exp3",
    width = 8,
    height = 5.5,
    caption_text = paste(
      "Figure 4. Brake-wear size comparison under equal mass loading.",
      "Model-predicted motile fraction as a function of log10(brake-wear concentration + 1) for coarse and fine brake-wear fractions.",
      "Points show raw observations and lines show model predictions."
    )
  )
  
  cat("Figure 4 saved.\n")
}

# ---------------------------------------------------------
# 12. Figure 5: Experiment 4
# ---------------------------------------------------------

exp4_raw <- dat |>
  dplyr::filter(experiment_num == 8.2)

if (!is.null(pred_exp4) && nrow(exp4_raw) > 0) {
  pred_exp4 <- pred_exp4 |>
    dplyr::mutate(
      days_from_start = suppressWarnings(as.integer(days_from_start)),
      day_label = paste0("Day ", days_from_start)
    )
  
  exp4_raw <- exp4_raw |>
    dplyr::mutate(
      day_label = paste0("Day ", days_from_start)
    )
  
  p_fig5 <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = exp4_raw,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = motility_ratio
      ),
      alpha = 0.30,
      size = 1.8,
      colour = particle_palette_local["SPM"]
    ) +
    ggplot2::geom_line(
      data = pred_exp4,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = fit_response
      ),
      linewidth = 0.9,
      colour = particle_palette_local["SPM"]
    ) +
    ggplot2::facet_wrap(~ day_label) +
    ggplot2::labs(
      title = "Field-derived SPM gradient",
      subtitle = "Model-predicted motile fraction versus log10(NTU + 1)",
      x = "log10(NTU + 1)",
      y = "Motile fraction",
      caption = "Points show raw observations and lines show model predictions."
    )
  
  manifest_figures[[length(manifest_figures) + 1]] <- save_plot_multi(
    plot = p_fig5,
    figure_name = "Fig5_field_spm_model_predictions",
    subdir = "models_exp4",
    width = 8.5,
    height = 6.5,
    caption_text = paste(
      "Figure 5. Field-derived SPM gradient.",
      "Model-predicted motile fraction as a function of log10(NTU + 1) across the full dataset (Days 4-12), facetted by day.",
      "Points show raw observations and lines show model predictions."
    )
  )
  
  cat("Figure 5 saved.\n")
}

# ---------------------------------------------------------
# 13. Supplementary Day 4 sensitivity figure
# ---------------------------------------------------------

if (!is.null(pred_exp4_day4) && nrow(exp4_raw) > 0) {
  exp4_day4_raw <- exp4_raw |>
    dplyr::filter(days_from_start == 4)
  
  p_figsx <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = exp4_day4_raw,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = motility_ratio
      ),
      alpha = 0.30,
      size = 1.8,
      colour = particle_palette_local["SPM"]
    ) +
    ggplot2::geom_line(
      data = pred_exp4_day4,
      ggplot2::aes(
        x = log10_ntu_plus1,
        y = fit_response
      ),
      linewidth = 0.9,
      colour = particle_palette_local["SPM"]
    ) +
    ggplot2::labs(
      title = "Field-derived SPM gradient: Day 4 sensitivity analysis",
      subtitle = "Model-predicted motile fraction versus log10(NTU + 1)",
      x = "log10(NTU + 1)",
      y = "Motile fraction",
      caption = "Points show raw observations and lines show model predictions."
    )
  
  manifest_figures[[length(manifest_figures) + 1]] <- save_plot_multi(
    plot = p_figsx,
    figure_name = "FigSx_field_spm_day4_sensitivity",
    subdir = "models_exp4",
    width = 8,
    height = 5.5,
    caption_text = paste(
      "Figure Sx. Experiment 4 Day 4 sensitivity analysis.",
      "Model-predicted motile fraction as a function of log10(NTU + 1) for the Day 4 subset only.",
      "Points show raw observations and lines show model predictions."
    )
  )
  
  cat("Supplementary Day 4 sensitivity figure saved.\n")
}

# ---------------------------------------------------------
# 14. Save manifest and log
# ---------------------------------------------------------

manifest_figures_df <- dplyr::bind_rows(manifest_figures)

file_manifest <- file.path(project_root, "outputs", "figures", "10_figures_main_manifest.csv")
readr::write_csv(manifest_figures_df, file_manifest)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_log <- file.path(dir_logs, paste0("10_figures_main_log_", timestamp_now, ".txt"))

sink(file_log)
cat("SPM Analysis - 10_figures_main log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")
cat("Input dataset:\n")
cat(input_source_used, "\n\n")
cat("Figures saved:\n")
print(manifest_figures_df)
cat("\n")
cat("Session information:\n\n")
print(utils::sessionInfo())
sink()

cat("Figure manifest written to:\n")
cat(file_manifest, "\n")
cat("Log written to:\n")
cat(file_log, "\n\n")

cat("========================================================\n")
cat("SCRIPT 10 COMPLETE: FIGURES MAIN\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")