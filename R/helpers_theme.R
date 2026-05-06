# =========================================================
# Script title: helpers_theme.R
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Project: SPM Analysis
# Manuscript: Manuscript 4
# Purpose: Define the central project plotting theme,
#          colour palettes, figure dimensions, and helper
#          functions for consistent publication-ready plots.
# Inputs: None
# Outputs: Global helper objects and functions available
#          to downstream scripts after sourcing
# Date created: 24 March 2026
# Last updated: 24 March 2026
# Notes/dependencies:
# - Source via 01_setup_packages_and_paths.R
# - Species colours are retained for cross-manuscript
#   consistency, but Article 4 figures should use
#   treatment-based colours rather than species colours.
# - Hard-coded hex values should not be used in plotting
#   scripts; access colours through named palette objects.
# =========================================================

message("Loading helper script: R/helpers_theme.R")

# ---------------------------------------------------------
# 1. Package checks
# ---------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required by helpers_theme.R.", call. = FALSE)
}

# ---------------------------------------------------------
# 2. Master colour palettes
# ---------------------------------------------------------

# Cross-manuscript species palette (fixed reference colours)
kelp_species_palette <- c(
  saccharina_latissima = "#ba89c3",  # mauve
  laminaria_digitata   = "#c09c0e",  # golden ochre
  laminaria_hyperborea = "#687fa2"   # steel blue
)

# Article 4 treatment palette
# Use these for plotting Manuscript 4 figures.
kelp_treatment_palette <- c(
  LIGHT_LOW   = "#a8c686",
  LIGHT_HIGH  = "#2f7d32",
  SAND        = "#f4e3b2",
  KAOLINITE   = "#d9d9d9",
  PEAT        = "#6b4226",
  BWC         = "#a64942",
  BWF         = "#5b2c2c",
  SPM         = "#6c6c6c",
  CONTROL     = "#4d4d4d"
)

# Extended project palette for flexible use in annotations,
# backgrounds, and secondary styling.
kelp_palette <- c(
  # species anchors
  SACCHARINA_LATISSIMA = kelp_species_palette[["saccharina_latissima"]],
  LAMINARIA_DIGITATA   = kelp_species_palette[["laminaria_digitata"]],
  LAMINARIA_HYPERBOREA = kelp_species_palette[["laminaria_hyperborea"]],

  # treatment keys
  LIGHT_LOW   = kelp_treatment_palette[["LIGHT_LOW"]],
  LIGHT_HIGH  = kelp_treatment_palette[["LIGHT_HIGH"]],
  SAND        = kelp_treatment_palette[["SAND"]],
  KAOLINITE   = kelp_treatment_palette[["KAOLINITE"]],
  PEAT        = kelp_treatment_palette[["PEAT"]],
  BWC         = kelp_treatment_palette[["BWC"]],
  BWF         = kelp_treatment_palette[["BWF"]],
  SPM         = kelp_treatment_palette[["SPM"]],
  CONTROL     = kelp_treatment_palette[["CONTROL"]],

  # neutrals / support
  BACKGROUND  = "#ffffff",
  GRID        = "#d9d9d9",
  AXIS        = "#2b2b2b",
  TEXT        = "#1f1f1f"
)

# ---------------------------------------------------------
# 3. Preferred mappings for common variables
# ---------------------------------------------------------

# Particle / treatment mapping for Article 4
kelp_particle_values <- c(
  "SAND"      = kelp_palette[["SAND"]],
  "KAOLINITE" = kelp_palette[["KAOLINITE"]],
  "PEAT"      = kelp_palette[["PEAT"]],
  "BWC"       = kelp_palette[["BWC"]],
  "BWF"       = kelp_palette[["BWF"]],
  "SPM"       = kelp_palette[["SPM"]],
  "CONTROL"   = kelp_palette[["CONTROL"]]
)

# Light gradient mapping for Experiment 1 grouped displays
kelp_light_values <- c(
  "LIGHT_LOW"  = kelp_palette[["LIGHT_LOW"]],
  "LIGHT_HIGH" = kelp_palette[["LIGHT_HIGH"]],
  "CONTROL"    = kelp_palette[["CONTROL"]]
)

# Optional linetype and shape mappings to support accessibility
kelp_linetype_values <- c(
  "SAND"      = "solid",
  "KAOLINITE" = "dashed",
  "PEAT"      = "dotted",
  "BWC"       = "solid",
  "BWF"       = "dashed",
  "SPM"       = "dotdash",
  "CONTROL"   = "solid"
)

kelp_shape_values <- c(
  "SAND"      = 21,
  "KAOLINITE" = 22,
  "PEAT"      = 24,
  "BWC"       = 21,
  "BWF"       = 22,
  "SPM"       = 24,
  "CONTROL"   = 19
)
# ---------------------------------------------------------
# 3b. Article 4 Day 4 figure colour families
# ---------------------------------------------------------

# Experiment-level anchor colours for synthesis figures
# such as Figure 8. Each experiment has one visual identity.
kelp_experiment_values <- c(
  "Experiment 1: light-only"        = kelp_palette[["CONTROL"]],
  "Experiment 2: defined particles" = kelp_palette[["LAMINARIA_DIGITATA"]],
  "Experiment 3: brake-wear size"   = kelp_palette[["SACCHARINA_LATISSIMA"]],
  "Experiment 4: field-derived SPM" = kelp_palette[["LIGHT_LOW"]]
)

# Experiment 1: light-only family.
# Kept deliberately neutral because the x-axis already carries
# the light-treatment structure.
kelp_exp1_light_family <- c(
  "0 lux"   = kelp_palette[["CONTROL"]],
  "4 lux"   = kelp_palette[["SPM"]],
  "70 lux"  = "#8fae78",
  "117 lux" = "#2f7d32"
)

# Experiment 2: defined-particle family.
# Earth/sediment tones are used for the defined-particle series.
kelp_exp2_particle_family <- c(
  "Sand"      = "#d8bd76",
  "Kaolinite" = kelp_palette[["LAMINARIA_DIGITATA"]],
  "Peat"      = kelp_palette[["PEAT"]]
)

# Experiment 3: brake-wear family.
# Coarse and fine fractions are kept in the same impact family,
# with the fine fraction shown as the darker/stronger colour.
kelp_exp3_brake_family <- c(
  "Coarse" = kelp_palette[["SACCHARINA_LATISSIMA"]],
  "Fine"   = kelp_palette[["BWF"]]
)

# Experiment 4: field-SPM family.
# Field-derived SPM is shown in the blue/green environmental family.
kelp_exp4_spm_family <- c(
  "SPM" = kelp_palette[["LIGHT_LOW"]]
)

# Shared support colours for model-output plots
kelp_raw_point_colour <- "#b8b8b8"
kelp_reference_line_colour <- kelp_palette[["CONTROL"]]
# ---------------------------------------------------------
# 4. Figure dimension constants
# ---------------------------------------------------------

# Journal-style widths stored in inches
fig_width_one_col <- 85 / 25.4
fig_width_two_col <- 178 / 25.4
fig_height_std    <- 110 / 25.4
fig_height_tall   <- 140 / 25.4
fig_height_short  <- 90 / 25.4
fig_dpi           <- 600

# Common point and line sizes
kelp_base_size        <- 11
kelp_linewidth        <- 0.5
kelp_panel_linewidth  <- 0.4
kelp_point_size       <- 2.2
kelp_raw_alpha        <- 0.35
kelp_ci_alpha         <- 0.18

# ---------------------------------------------------------
# 5. Core theme object
# ---------------------------------------------------------

theme_kelp <- function(base_size = kelp_base_size,
                       base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1,
        colour = kelp_palette[["TEXT"]],
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size,
        colour = kelp_palette[["TEXT"]],
        hjust = 0
      ),
      plot.caption = ggplot2::element_text(
        size = base_size - 2,
        colour = kelp_palette[["TEXT"]],
        hjust = 0
      ),
      axis.title = ggplot2::element_text(
        size = base_size,
        colour = kelp_palette[["TEXT"]]
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1,
        colour = kelp_palette[["AXIS"]]
      ),
      axis.line = ggplot2::element_line(
        colour = kelp_palette[["AXIS"]],
        linewidth = kelp_panel_linewidth
      ),
      axis.ticks = ggplot2::element_line(
        colour = kelp_palette[["AXIS"]],
        linewidth = kelp_panel_linewidth
      ),
      panel.grid.major = ggplot2::element_line(
        colour = kelp_palette[["GRID"]],
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = kelp_palette[["AXIS"]],
        fill = NA,
        linewidth = kelp_panel_linewidth
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey92",
        colour = "grey70",
        linewidth = kelp_panel_linewidth
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = kelp_palette[["TEXT"]],
        size = base_size - 1
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        face = "bold",
        size = base_size - 1,
        colour = kelp_palette[["TEXT"]]
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 1,
        colour = kelp_palette[["TEXT"]]
      ),
      legend.key = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = kelp_palette[["BACKGROUND"]],
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = kelp_palette[["BACKGROUND"]],
        colour = NA
      )
    )
}

# ---------------------------------------------------------
# 6. Convenience theme setter
# ---------------------------------------------------------

set_kelp_theme <- function(base_size = kelp_base_size,
                           base_family = "") {
  ggplot2::theme_set(theme_kelp(base_size = base_size, base_family = base_family))
  invisible(TRUE)
}

# ---------------------------------------------------------
# 7. Validation helpers
# ---------------------------------------------------------

validate_palette_keys <- function(keys, palette) {
  missing_keys <- setdiff(keys, names(palette))
  if (length(missing_keys) > 0) {
    stop(
      paste0(
        "Unknown palette key(s): ",
        paste(missing_keys, collapse = ", "),
        ". Valid keys are: ",
        paste(names(palette), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

get_kelp_cols <- function(keys, palette = kelp_palette) {
  validate_palette_keys(keys, palette)
  unname(palette[keys])
}

# ---------------------------------------------------------
# 8. Manual scale helpers
# ---------------------------------------------------------

scale_colour_kelp_experiment <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_experiment_values,
    ...
  )
}

scale_fill_kelp_experiment <- function(...) {
  ggplot2::scale_fill_manual(
    values = kelp_experiment_values,
    ...
  )
}

scale_colour_kelp_exp1_light <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_exp1_light_family,
    ...
  )
}

scale_colour_kelp_exp2_particle <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_exp2_particle_family,
    ...
  )
}

scale_colour_kelp_exp3_brake <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_exp3_brake_family,
    ...
  )
}

scale_fill_kelp_exp4_spm <- function(...) {
  ggplot2::scale_fill_manual(
    values = kelp_exp4_spm_family,
    ...
  )
}

scale_colour_kelp_exp4_spm <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_exp4_spm_family,
    ...
  )
}

scale_colour_kelp_particle <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_particle_values,
    ...
  )
}

scale_fill_kelp_particle <- function(...) {
  ggplot2::scale_fill_manual(
    values = kelp_particle_values,
    ...
  )
}

scale_colour_kelp_light <- function(...) {
  ggplot2::scale_colour_manual(
    values = kelp_light_values,
    ...
  )
}

scale_fill_kelp_light <- function(...) {
  ggplot2::scale_fill_manual(
    values = kelp_light_values,
    ...
  )
}

scale_linetype_kelp <- function(...) {
  ggplot2::scale_linetype_manual(
    values = kelp_linetype_values,
    ...
  )
}

scale_shape_kelp <- function(...) {
  ggplot2::scale_shape_manual(
    values = kelp_shape_values,
    ...
  )
}

# Flexible manual scales using any named subset of kelp_palette
scale_colour_kelp_manual <- function(values, ...) {
  validate_palette_keys(values, kelp_palette)
  ggplot2::scale_colour_manual(
    values = kelp_palette[values],
    ...
  )
}

scale_fill_kelp_manual <- function(values, ...) {
  validate_palette_keys(values, kelp_palette)
  ggplot2::scale_fill_manual(
    values = kelp_palette[values],
    ...
  )
}

# ---------------------------------------------------------
# 9. Label helpers for common treatment names
# ---------------------------------------------------------

label_particle_type <- function(x) {
  dplyr::recode(
    x,
    "SAND"      = "Sand",
    "KAOLINITE" = "Kaolinite",
    "PEAT"      = "Peat",
    "BWC"       = "Brake wear coarse",
    "BWF"       = "Brake wear fine",
    "SPM"       = "Field SPM",
    "CONTROL"   = "Control",
    .default = x
  )
}

label_size_class <- function(x) {
  dplyr::recode(
    x,
    "BWC"    = "Coarse",
    "BWF"    = "Fine",
    "COARSE" = "Coarse",
    "FINE"   = "Fine",
    .default = x
  )
}

# ---------------------------------------------------------
# 10. Plot finishing helper
# ---------------------------------------------------------

finish_kelp_plot <- function(p,
                             title = NULL,
                             subtitle = NULL,
                             caption = NULL) {
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    theme_kelp()
}

# ---------------------------------------------------------
# 11. Recommended default theme on source
# ---------------------------------------------------------

set_kelp_theme()

message("helpers_theme.R loaded successfully.")
message("Available palette objects: kelp_species_palette, kelp_treatment_palette, kelp_palette")
message("Available theme helpers: theme_kelp(), set_kelp_theme(), finish_kelp_plot()")
message("Available scale helpers: scale_colour_kelp_particle(), scale_fill_kelp_particle(), scale_colour_kelp_light(), scale_fill_kelp_light()")
message("Available Day 4 palettes: kelp_experiment_values, kelp_exp1_light_family, kelp_exp2_particle_family, kelp_exp3_brake_family, kelp_exp4_spm_family")
message("Available Day 4 scale helpers: scale_colour_kelp_experiment(), scale_colour_kelp_exp1_light(), scale_colour_kelp_exp2_particle(), scale_colour_kelp_exp3_brake(), scale_colour_kelp_exp4_spm()")