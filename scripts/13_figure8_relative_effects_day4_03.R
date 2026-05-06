# =========================================================
# Script title: 13_figure8_relative_effects_day4_02.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4 / Article 4
# Purpose: Assemble Day 4 relative-effect outputs from
#          Experiments 1-4 and generate Figure 8, a synthesis
#          forest plot comparing predicted proportional changes
#          in zoospore motility across the experimental series.
#
# Inputs:
# - outputs/tables/day4/06b_exp1_light_day4_relative_effects.csv
# - outputs/tables/day4/07b_exp2_defined_particles_day4_relative_effects.csv
# - outputs/tables/day4/08b_exp3_brake_size_day4_relative_effects.csv
# - outputs/tables/day4/09b_exp4_field_spm_day4_relative_effects.csv
#
# Outputs:
# - outputs/tables/day4/13_figure8_day4_relative_effects_combined.csv
# - outputs/tables/day4/13_figure8_day4_effect_ranking.csv
# - outputs/tables/day4/13_figure8_day4_reference_conditions.csv
# - outputs/figures/day4/Fig8_day4_relative_effects_forest_plot.{pdf,png,tiff}
# - outputs/logs/day4/13_figure8_relative_effects_day4_log_*.txt
#
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 06b, 07b, 08b, and 09b Day 4 model scripts before
#   this script.
# - This script does not fit models.
# - Effects are predicted proportional changes in motile
#   fraction relative to experiment-specific reference
#   conditions.
# - Reference conditions:
#   Experiment 1: 117 lux
#   Experiment 2: Sand at 25 NTU
#   Experiment 3: BWC at 4.5 ug/L
#   Experiment 4: lowest observed NTU
# - Experiment 2 is collapsed to one row per particle type
#   because the preferred model was particle-only and repeated
#   NTU-level contrasts are duplicates for Figure 8.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 13: FIGURE 8 DAY 4 RELATIVE EFFECTS\n")
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
  "tibble",
  "forcats"
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

dir_day4_tables <- file.path(project_root, "outputs", "tables", "day4")
dir_day4_logs <- file.path(project_root, "outputs", "logs", "day4")
dir_day4_figures <- file.path(project_root, "outputs", "figures", "day4")

dir.create(dir_day4_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_day4_figures, recursive = TRUE, showWarnings = FALSE)

file_exp1_effects <- file.path(
  dir_day4_tables,
  "06b_exp1_light_day4_relative_effects.csv"
)

file_exp2_effects <- file.path(
  dir_day4_tables,
  "07b_exp2_defined_particles_day4_relative_effects.csv"
)

file_exp3_effects <- file.path(
  dir_day4_tables,
  "08b_exp3_brake_size_day4_relative_effects.csv"
)

file_exp4_effects <- file.path(
  dir_day4_tables,
  "09b_exp4_field_spm_day4_relative_effects.csv"
)

file_combined <- file.path(
  dir_day4_tables,
  "13_figure8_day4_relative_effects_combined.csv"
)

file_ranking <- file.path(
  dir_day4_tables,
  "13_figure8_day4_effect_ranking.csv"
)

file_references <- file.path(
  dir_day4_tables,
  "13_figure8_day4_reference_conditions.csv"
)

file_manifest <- file.path(
  dir_day4_tables,
  "13_figure8_day4_manifest.csv"
)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")

file_log <- file.path(
  dir_day4_logs,
  paste0("13_figure8_relative_effects_day4_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 5. Small local helper functions
# ---------------------------------------------------------

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      paste0(
        "Required input file not found for ",
        label,
        ":\n",
        path
      ),
      call. = FALSE
    )
  }
  
  readr::read_csv(path, show_col_types = FALSE)
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
                                  width = 9,
                                  height = 8,
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

check_required_effect_columns <- function(df, label) {
  required_cols <- c(
    "experiment",
    "experiment_num",
    "experiment_label",
    "analysis_branch",
    "model",
    "response",
    "contrast_label",
    "reference_condition",
    "treatment_condition",
    "fit_control",
    "fit_treatment",
    "percent_change",
    "conf_low_percent",
    "conf_high_percent",
    "absolute_change"
  )
  
  missing_cols <- required_cols[!required_cols %in% names(df)]
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        label,
        " relative-effects file is missing required column(s):\n- ",
        paste(missing_cols, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

# ---------------------------------------------------------
# 6. Load relative-effect outputs
# ---------------------------------------------------------

exp1_effects <- read_required_csv(file_exp1_effects, "Experiment 1")
exp2_effects <- read_required_csv(file_exp2_effects, "Experiment 2")
exp3_effects <- read_required_csv(file_exp3_effects, "Experiment 3")
exp4_effects <- read_required_csv(file_exp4_effects, "Experiment 4")

check_required_effect_columns(exp1_effects, "Experiment 1")
check_required_effect_columns(exp2_effects, "Experiment 2")
check_required_effect_columns(exp3_effects, "Experiment 3")
check_required_effect_columns(exp4_effects, "Experiment 4")

cat("Relative-effect files loaded and verified.\n")
cat("Experiment 1 rows:", nrow(exp1_effects), "\n")
cat("Experiment 2 rows:", nrow(exp2_effects), "\n")
cat("Experiment 3 rows:", nrow(exp3_effects), "\n")
cat("Experiment 4 rows:", nrow(exp4_effects), "\n\n")

# ---------------------------------------------------------
# 7. Harmonise Experiment 1 labels
# ---------------------------------------------------------

exp1_fig8 <- exp1_effects |>
  dplyr::mutate(
    experiment_block = "Experiment 1: light-only",
    experiment_order = 1L,
    display_label = dplyr::case_when(
      stringr::str_detect(contrast_label, "^0 lux") ~ "0 lux",
      stringr::str_detect(contrast_label, "^4 lux") ~ "4 lux",
      stringr::str_detect(contrast_label, "^70 lux") ~ "70 lux",
      TRUE ~ treatment_condition
    ),
    contrast_group = "Light treatment",
    reference_short = "117 lux",
    include_in_fig8 = TRUE,
    fig8_note = "All Experiment 1 contrasts retained."
  )

# ---------------------------------------------------------
# 8. Harmonise and collapse Experiment 2 labels
# ---------------------------------------------------------

# The Experiment 2 preferred model was particle-only, so
# the NTU-specific kaolinite and peat rows are duplicates for
# Figure 8. Collapse to one representative row per particle
# type: Kaolinite and Peat relative to Sand.
exp2_fig8 <- exp2_effects |>
  dplyr::mutate(
    particle_type_for_collapse = dplyr::case_when(
      "particle_type_plot" %in% names(exp2_effects) ~ as.character(particle_type_plot),
      stringr::str_detect(treatment_condition, regex("kaolinite", ignore_case = TRUE)) ~ "Kaolinite",
      stringr::str_detect(treatment_condition, regex("peat", ignore_case = TRUE)) ~ "Peat",
      stringr::str_detect(treatment_condition, regex("sand", ignore_case = TRUE)) ~ "Sand",
      TRUE ~ treatment_condition
    )
  ) |>
  dplyr::filter(
    particle_type_for_collapse %in% c("Kaolinite", "Peat")
  ) |>
  dplyr::group_by(particle_type_for_collapse) |>
  dplyr::arrange(abs(percent_change), .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    experiment_block = "Experiment 2: defined particles",
    experiment_order = 2L,
    display_label = dplyr::case_when(
      particle_type_for_collapse == "Kaolinite" ~ "Kaolinite",
      particle_type_for_collapse == "Peat" ~ "Peat",
      TRUE ~ treatment_condition
    ),
    contrast_group = "Particle type",
    reference_short = "Sand at 25 NTU",
    include_in_fig8 = TRUE,
    fig8_note = "Experiment 2 collapsed to one contrast per particle type because preferred model was particle-only."
  )

# ---------------------------------------------------------
# 9. Harmonise Experiment 3 labels
# ---------------------------------------------------------

exp3_fig8 <- exp3_effects |>
  dplyr::mutate(
    experiment_block = "Experiment 3: brake-wear size",
    experiment_order = 3L,
    display_label = dplyr::case_when(
      stringr::str_detect(treatment_condition, regex("Coarse at 450", ignore_case = TRUE)) ~ "Coarse 450 ug/L",
      stringr::str_detect(treatment_condition, regex("Coarse at 45", ignore_case = TRUE)) ~ "Coarse 45 ug/L",
      stringr::str_detect(treatment_condition, regex("Fine at 450", ignore_case = TRUE)) ~ "Fine 450 ug/L",
      stringr::str_detect(treatment_condition, regex("Fine at 45", ignore_case = TRUE)) ~ "Fine 45 ug/L",
      stringr::str_detect(treatment_condition, regex("Fine at 4.5", ignore_case = TRUE)) ~ "Fine 4.5 ug/L",
      TRUE ~ treatment_condition
    ),
    contrast_group = "Brake-wear fraction/load",
    reference_short = "BWC 4.5 ug/L",
    include_in_fig8 = TRUE,
    fig8_note = "All Experiment 3 contrasts retained because preferred model included size by load interaction."
  )

# ---------------------------------------------------------
# 10. Harmonise Experiment 4 labels
# ---------------------------------------------------------

exp4_fig8 <- exp4_effects |>
  dplyr::mutate(
    experiment_block = "Experiment 4: field-derived SPM",
    experiment_order = 4L,
    display_label = dplyr::case_when(
      "ntu" %in% names(exp4_effects) & ntu == 4144 ~ "4144 NTU (extreme)",
      "ntu" %in% names(exp4_effects) ~ paste0(round(ntu, 0), " NTU"),
      TRUE ~ treatment_condition
    ),
    contrast_group = "Field-SPM NTU",
    reference_short = "0 NTU",
    include_in_fig8 = TRUE,
    fig8_note = "Experiment 4 contrasts shown relative to lowest observed NTU."
  )
# ---------------------------------------------------------
# 11. Combine Figure 8 effects
# ---------------------------------------------------------

combined_effects <- dplyr::bind_rows(
  exp1_fig8,
  exp2_fig8,
  exp3_fig8,
  exp4_fig8
) |>
  dplyr::mutate(
    percent_change = as.numeric(percent_change),
    conf_low_percent = as.numeric(conf_low_percent),
    conf_high_percent = as.numeric(conf_high_percent),
    absolute_percent_change = abs(percent_change),

    effect_direction = dplyr::case_when(
      percent_change < 0 ~ "Decrease",
      percent_change > 0 ~ "Increase",
      TRUE ~ "No change"
    ),

    ci_crosses_zero = conf_low_percent <= 0 & conf_high_percent >= 0,

    effect_strength_label = dplyr::case_when(
      absolute_percent_change >= 25 ~ ">=25%",
      absolute_percent_change >= 10 ~ "10-25%",
      absolute_percent_change > 0 ~ "<10%",
      TRUE ~ "0%"
    ),

    display_label = as.character(display_label),

    experiment_block = factor(
      experiment_block,
      levels = c(
        "Experiment 1: light-only",
        "Experiment 2: defined particles",
        "Experiment 3: brake-wear size",
        "Experiment 4: field-derived SPM"
      )
    ),

    within_experiment_order = dplyr::case_when(
      # Experiment 1: light treatments
      experiment_block == "Experiment 1: light-only" &
        display_label == "0 lux" ~ 1,

      experiment_block == "Experiment 1: light-only" &
        display_label == "4 lux" ~ 2,

      experiment_block == "Experiment 1: light-only" &
        display_label == "70 lux" ~ 3,

      # Experiment 2: particle identity
      experiment_block == "Experiment 2: defined particles" &
        display_label == "Kaolinite" ~ 1,

      experiment_block == "Experiment 2: defined particles" &
        display_label == "Peat" ~ 2,

      # Experiment 3: brake-wear fraction/load
      experiment_block == "Experiment 3: brake-wear size" &
        display_label == "Coarse 45 ug/L" ~ 1,

      experiment_block == "Experiment 3: brake-wear size" &
        display_label == "Coarse 450 ug/L" ~ 2,

      experiment_block == "Experiment 3: brake-wear size" &
        display_label == "Fine 4.5 ug/L" ~ 3,

      experiment_block == "Experiment 3: brake-wear size" &
        display_label == "Fine 45 ug/L" ~ 4,

      experiment_block == "Experiment 3: brake-wear size" &
        display_label == "Fine 450 ug/L" ~ 5,

      # Experiment 4: field-SPM NTU gradient
      experiment_block == "Experiment 4: field-derived SPM" &
        display_label == "25 NTU" ~ 1,

      experiment_block == "Experiment 4: field-derived SPM" &
        display_label == "100 NTU" ~ 2,

      experiment_block == "Experiment 4: field-derived SPM" &
        display_label == "400 NTU" ~ 3,

      experiment_block == "Experiment 4: field-derived SPM" &
        display_label == "500 NTU" ~ 4,

      experiment_block == "Experiment 4: field-derived SPM" &
        display_label == "4144 NTU (extreme)" ~ 5,

      TRUE ~ 999
    )
  ) |>
  dplyr::arrange(
    experiment_order,
    within_experiment_order
  ) |>
  dplyr::group_by(experiment_block) |>
  dplyr::mutate(
    display_order = dplyr::row_number()
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    display_label_unique = paste0(
      display_label,
      "___",
      as.integer(experiment_order),
      "_",
      sprintf("%02d", display_order)
    ),

    # ggplot places the first factor level at the bottom,
    # so reverse the ordered labels to display the intended
    # sequence from top to bottom within each facet.
    display_label_unique = factor(
      display_label_unique,
      levels = rev(unique(display_label_unique))
    )
  )

readr::write_csv(combined_effects, file_combined)

cat("Combined Figure 8 relative effects written to:\n")
cat(file_combined, "\n\n")

print(
  combined_effects |>
    dplyr::select(
      experiment_block,
      display_label,
      within_experiment_order,
      reference_short,
      percent_change,
      conf_low_percent,
      conf_high_percent,
      ci_crosses_zero,
      model
    )
)

# ---------------------------------------------------------
# 12. Reference condition summary
# ---------------------------------------------------------

reference_conditions <- combined_effects |>
  dplyr::distinct(
    experiment,
    experiment_num,
    experiment_label,
    experiment_block,
    reference_condition,
    reference_short
  ) |>
  dplyr::arrange(experiment_num)

readr::write_csv(reference_conditions, file_references)

cat("Reference condition summary written to:\n")
cat(file_references, "\n\n")

print(reference_conditions)

# ---------------------------------------------------------
# 13. Effect ranking summary
# ---------------------------------------------------------

largest_negative <- combined_effects |>
  dplyr::filter(percent_change == min(percent_change, na.rm = TRUE)) |>
  dplyr::slice(1) |>
  dplyr::mutate(ranking_category = "largest_negative_effect")

largest_positive <- combined_effects |>
  dplyr::filter(percent_change == max(percent_change, na.rm = TRUE)) |>
  dplyr::slice(1) |>
  dplyr::mutate(ranking_category = "largest_positive_effect")

largest_absolute <- combined_effects |>
  dplyr::filter(absolute_percent_change == max(absolute_percent_change, na.rm = TRUE)) |>
  dplyr::slice(1) |>
  dplyr::mutate(ranking_category = "largest_absolute_effect")

effect_ranking <- dplyr::bind_rows(
  largest_negative,
  largest_positive,
  largest_absolute
) |>
  dplyr::select(
    ranking_category,
    experiment,
    experiment_num,
    experiment_label,
    experiment_block,
    display_label,
    contrast_label,
    reference_condition,
    treatment_condition,
    model,
    fit_control,
    fit_treatment,
    percent_change,
    conf_low_percent,
    conf_high_percent,
    absolute_percent_change,
    ci_crosses_zero
  )

readr::write_csv(effect_ranking, file_ranking)

cat("Effect ranking summary written to:\n")
cat(file_ranking, "\n\n")

print(effect_ranking)

# ---------------------------------------------------------
# 14. Figure 8 forest plot
# ---------------------------------------------------------

if (exists("set_kelp_theme", mode = "function", inherits = TRUE)) {
  set_kelp_theme()
}

experiment_palette <- c(
  "Experiment 1: light-only" = "#6E6E6E",
  "Experiment 2: defined particles" = "#A6761D",
  "Experiment 3: brake-wear size" = "#7570B3",
  "Experiment 4: field-derived SPM" = "#2C7F62"
)

# Clean axis labels by removing unique suffix.
clean_y_labels <- function(x) {
  stringr::str_replace(as.character(x), "___.*$", "")
}

# X-axis range based on observed effects, with padding.
x_min <- min(combined_effects$conf_low_percent, na.rm = TRUE)
x_max <- max(combined_effects$conf_high_percent, na.rm = TRUE)

x_pad <- max(5, 0.08 * (x_max - x_min))

x_limits <- c(
  floor((x_min - x_pad) / 10) * 10,
  ceiling((x_max + x_pad) / 10) * 10
)

p_fig8 <- ggplot2::ggplot(
  combined_effects,
  ggplot2::aes(
    x = percent_change,
    y = display_label_unique
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    colour = "grey35"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = conf_low_percent,
      xmax = conf_high_percent,
      colour = experiment_block
    ),
    height = 0.18,
    linewidth = 0.7,
    alpha = 0.9
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = experiment_block,
      shape = ci_crosses_zero
    ),
    size = 2.8
  ) +
  ggplot2::facet_grid(
    experiment_block ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_y_discrete(labels = clean_y_labels) +
  ggplot2::scale_colour_manual(values = experiment_palette, drop = FALSE) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 16),
    labels = c(`TRUE` = "CI crosses 0", `FALSE` = "CI excludes 0"),
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(xlim = x_limits) +
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
      "Experiment 2 is collapsed to one contrast per particle type because the preferred model was particle-only."
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
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8.5
    )
  )

figure_manifest <- save_plot_multi_local(
  plot = p_fig8,
  file_stem = "Fig8_day4_relative_effects_forest_plot",
  width = 9.5,
  height = 8.5,
  dpi = 600
)

cat("Figure 8 forest plot saved:\n")
print(figure_manifest)

# ---------------------------------------------------------
# 15. Write manifest
# ---------------------------------------------------------

manifest <- tibble::tibble(
  output_type = c(
    "combined_effects",
    "effect_ranking",
    "reference_conditions",
    "figure8_forest_plot"
  ),
  path = c(
    file_combined,
    file_ranking,
    file_references,
    figure_manifest$file_pdf[1]
  ),
  notes = c(
    "Combined harmonised relative effects from Experiments 1-4.",
    "Largest negative, largest positive, and largest absolute effects.",
    "Experiment-specific reference conditions used for proportional effects.",
    "PDF version of Figure 8 forest plot; PNG and TIFF also saved."
  )
)

readr::write_csv(manifest, file_manifest)

cat("Manifest written to:\n")
cat(file_manifest, "\n\n")

print(manifest)

# ---------------------------------------------------------
# 16. Write log
# ---------------------------------------------------------

sink(file_log)

cat("Figure 8 Day 4 relative effects log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input files:\n")
cat("Experiment 1:", file_exp1_effects, "\n")
cat("Experiment 2:", file_exp2_effects, "\n")
cat("Experiment 3:", file_exp3_effects, "\n")
cat("Experiment 4:", file_exp4_effects, "\n\n")

cat("Rows loaded:\n")
cat("Experiment 1:", nrow(exp1_effects), "\n")
cat("Experiment 2:", nrow(exp2_effects), "\n")
cat("Experiment 3:", nrow(exp3_effects), "\n")
cat("Experiment 4:", nrow(exp4_effects), "\n\n")

cat("Rows included in Figure 8 combined dataset:", nrow(combined_effects), "\n\n")

cat("Reference conditions:\n")
print(reference_conditions)

cat("\nEffect ranking:\n")
print(effect_ranking)

cat("\nCombined effects:\n")
print(
  combined_effects |>
    dplyr::select(
      experiment_block,
      display_label,
      percent_change,
      conf_low_percent,
      conf_high_percent,
      ci_crosses_zero,
      model
    )
)

cat("\nFigure manifest:\n")
print(figure_manifest)

sink()

cat("Log written to:\n")
cat(file_log, "\n\n")

cat("SCRIPT 13 COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")

