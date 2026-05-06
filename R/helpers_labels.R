# =========================================================
# Script title: helpers_labels.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Define central label helpers for experiments,
#          treatments, variables, QC outputs, figures,
#          tables, and publication-ready wording.
# Inputs: Character vectors, variable names, factor levels
# Outputs: Reusable helper functions and named label maps
#          for consistent wording across the project
# Author: Marianne Glascott
# Date created: 25 February 2026
# Last updated: 24 March 2026
# Notes/dependencies:
# - Source via 01_setup_packages_and_paths.R
# - Prevents duplicated relabelling logic
#   across analysis, figure, and table scripts.
# - Labels are manuscript-facing where appropriate.
# =========================================================

message("Loading helper script: R/helpers_labels.R")

# ---------------------------------------------------------
# 1. Package checks
# ---------------------------------------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required by helpers_labels.R.", call. = FALSE)
}

# ---------------------------------------------------------
# 2. Core label dictionaries
# ---------------------------------------------------------

experiment_labels <- c(
  "LIGHT_ONLY"        = "Experiment 1: Light-only gradient",
  "EXP1"              = "Experiment 1: Light-only gradient",
  "LIGHT"             = "Experiment 1: Light-only gradient",
  "DEFINED_PARTICLES" = "Experiment 2: Defined particle concentration series",
  "EXP2"              = "Experiment 2: Defined particle concentration series",
  "BRAKE_SIZE"        = "Experiment 3: Particle size comparison under equal mass loading",
  "EXP3"              = "Experiment 3: Particle size comparison under equal mass loading",
  "FIELD_SPM"         = "Experiment 4: Field-derived SPM gradient",
  "EXP4"              = "Experiment 4: Field-derived SPM gradient"
)

experiment_short_labels <- c(
  "LIGHT_ONLY"        = "Exp. 1 Light-only",
  "EXP1"              = "Exp. 1 Light-only",
  "LIGHT"             = "Exp. 1 Light-only",
  "DEFINED_PARTICLES" = "Exp. 2 Defined particles",
  "EXP2"              = "Exp. 2 Defined particles",
  "BRAKE_SIZE"        = "Exp. 3 Brake size",
  "EXP3"              = "Exp. 3 Brake size",
  "FIELD_SPM"         = "Exp. 4 Field SPM",
  "EXP4"              = "Exp. 4 Field SPM"
)

particle_labels <- c(
  "SAND"      = "Sand",
  "KAOLINITE" = "Kaolinite",
  "PEAT"      = "Peat",
  "BWC"       = "Brake wear coarse",
  "BWF"       = "Brake wear fine",
  "SPM"       = "Field SPM",
  "CONTROL"   = "Control"
)

particle_short_labels <- c(
  "SAND"      = "Sand",
  "KAOLINITE" = "Kaolinite",
  "PEAT"      = "Peat",
  "BWC"       = "BWC",
  "BWF"       = "BWF",
  "SPM"       = "SPM",
  "CONTROL"   = "Control"
)

size_class_labels <- c(
  "COARSE" = "Coarse",
  "FINE"   = "Fine",
  "BWC"    = "Coarse",
  "BWF"    = "Fine"
)

season_labels <- c(
  "SPRING" = "Spring",
  "SUMMER" = "Summer",
  "AUTUMN" = "Autumn",
  "FALL"   = "Autumn",
  "WINTER" = "Winter"
)

species_labels <- c(
  "Laminaria digitata"   = "Laminaria digitata",
  "Laminaria hyperborea" = "Laminaria hyperborea",
  "Saccharina latissima" = "Saccharina latissima",
  "LD" = "Laminaria digitata",
  "LH" = "Laminaria hyperborea",
  "SL" = "Saccharina latissima"
)

# ---------------------------------------------------------
# 3. Variable labels for axes, tables, and outputs
# ---------------------------------------------------------

variable_labels <- c(
  "mobile_cell_count"      = "Mobile cell count",
  "stationary_cell_count"  = "Stationary cell count",
  "total_cells"            = "Total cell count",
  "motility_ratio"         = "Motile fraction",
  "days_from_start"        = "Days from start",
  "lux_exposure"           = "Irradiance",
  "ntu"                    = "NTU",
  "mass_loading_ug_l"      = expression("Mass loading ("*mu*"g L"^-1*")"),
  "cu_ug_l"                = expression("Copper ("*mu*"g L"^-1*")"),
  "particle_type"          = "Particle type",
  "particle_class"         = "Particle class",
  "size_class"             = "Size class",
  "culture"                = "Culture",
  "well"                   = "Well",
  "experiment_id"          = "Experiment ID",
  "experiment_type"        = "Experiment",
  "count_date_yyyy_mm_dd"  = "Count date",
  "start_date_yyyy_mm_dd"  = "Start date",
  "frame_count"            = "Frame count",
  "tile_count"             = "Tile count"
)

variable_labels_plain <- c(
  "mobile_cell_count"      = "Mobile cell count",
  "stationary_cell_count"  = "Stationary cell count",
  "total_cells"            = "Total cell count",
  "motility_ratio"         = "Motile fraction",
  "days_from_start"        = "Days from start",
  "lux_exposure"           = "Irradiance",
  "ntu"                    = "NTU",
  "mass_loading_ug_l"      = "Mass loading (ug/L)",
  "cu_ug_l"                = "Copper (ug/L)",
  "particle_type"          = "Particle type",
  "particle_class"         = "Particle class",
  "size_class"             = "Size class",
  "culture"                = "Culture",
  "well"                   = "Well",
  "experiment_id"          = "Experiment ID",
  "experiment_type"        = "Experiment",
  "count_date_yyyy_mm_dd"  = "Count date",
  "start_date_yyyy_mm_dd"  = "Start date",
  "frame_count"            = "Frame count",
  "tile_count"             = "Tile count"
)

# ---------------------------------------------------------
# 4. QC flag and exclusion labels
# ---------------------------------------------------------

qc_flag_labels <- c(
  "flag_missing_response_counts"   = "Missing response counts",
  "flag_total_cell_count_le_zero"  = "Total cell count less than or equal to zero",
  "flag_inconsistent_count_totals" = "Inconsistent count totals",
  "flag_missing_key_predictor"     = "Missing key predictor for experiment",
  "flag_unexpected_frame_count"    = "Unexpected frame count",
  "flag_impossible_date_sequences" = "Impossible date sequence",
  "flag_negative_mobile_count"     = "Negative mobile count",
  "flag_negative_stationary_count" = "Negative stationary count",
  "flag_negative_ntu"              = "Negative NTU",
  "flag_negative_lux"              = "Negative irradiance",
  "flag_negative_mass_loading"     = "Negative mass loading"
)

exclusion_labels <- c(
  "exclude_missing_response_counts"   = "Excluded: missing response counts",
  "exclude_total_cell_count_le_zero"  = "Excluded: total cell count less than or equal to zero",
  "exclude_missing_key_predictor"     = "Excluded: missing key predictor for experiment",
  "exclude_negative_mobile_count"     = "Excluded: negative mobile count",
  "exclude_negative_stationary_count" = "Excluded: negative stationary count"
)

# ---------------------------------------------------------
# 5. Figure title helpers
# ---------------------------------------------------------

figure_title_labels <- c(
  "fig1_design"        = "Experimental design overview",
  "fig2_light"         = "Light-only gradient",
  "fig3_particles"     = "Defined particle concentration series",
  "fig4_brake_size"    = "Particle size comparison under equal mass loading",
  "fig5_field_spm"     = "Field-derived SPM gradient",
  "figS1_qc"           = "QC flow diagram or exclusion summary",
  "figS2_eda"          = "EDA distributions of response counts and motility ratio",
  "figS3_sensitivity"  = "Sensitivity model comparisons",
  "figS4_pilot"        = "Optional pilot mineral gradient"
)

# ---------------------------------------------------------
# 6. Table title helpers
# ---------------------------------------------------------

table_title_labels <- c(
  "table1_design"      = "Experimental design summary",
  "table2_dataset"     = "Dataset summary",
  "table3_models"      = "Model summary table",
  "table4_sensitivity" = "Model comparison / sensitivity summary"
)

# ---------------------------------------------------------
# 7. Safe label getter
# ---------------------------------------------------------

get_label <- function(x, dictionary, default_to_input = TRUE) {
  out <- unname(dictionary[as.character(x)])
  if (default_to_input) {
    out[is.na(out)] <- as.character(x)[is.na(out)]
  }
  out
}

# ---------------------------------------------------------
# 8. Vector relabelling helpers
# ---------------------------------------------------------

label_experiment <- function(x, short = FALSE) {
  dict <- if (short) experiment_short_labels else experiment_labels
  get_label(x, dict)
}

label_particle_type <- function(x, short = FALSE) {
  dict <- if (short) particle_short_labels else particle_labels
  get_label(x, dict)
}

label_size_class <- function(x) {
  get_label(x, size_class_labels)
}

label_season <- function(x) {
  get_label(x, season_labels)
}

label_species <- function(x) {
  get_label(x, species_labels)
}

label_variable <- function(x, plain = FALSE) {
  dict <- if (plain) variable_labels_plain else variable_labels
  out <- unname(dict[as.character(x)])
  missing_idx <- is.na(out)
  if (any(missing_idx)) {
    out[missing_idx] <- as.character(x)[missing_idx]
  }
  if (length(out) == 1) {
    return(out[[1]])
  }
  out
}

label_qc_flag <- function(x) {
  get_label(x, qc_flag_labels)
}

label_exclusion <- function(x) {
  get_label(x, exclusion_labels)
}

label_figure_title <- function(x) {
  get_label(x, figure_title_labels)
}

label_table_title <- function(x) {
  get_label(x, table_title_labels)
}

# ---------------------------------------------------------
# 9. Data frame relabelling helpers
# ---------------------------------------------------------

apply_experiment_labels <- function(data, column = "experiment_type", short = FALSE) {
  if (!column %in% names(data)) {
    stop(paste0("Column not found in data: ", column), call. = FALSE)
  }
  data[[column]] <- label_experiment(data[[column]], short = short)
  data
}

apply_particle_labels <- function(data, column = "particle_type", short = FALSE) {
  if (!column %in% names(data)) {
    stop(paste0("Column not found in data: ", column), call. = FALSE)
  }
  data[[column]] <- label_particle_type(data[[column]], short = short)
  data
}

apply_size_class_labels <- function(data, column = "size_class") {
  if (!column %in% names(data)) {
    stop(paste0("Column not found in data: ", column), call. = FALSE)
  }
  data[[column]] <- label_size_class(data[[column]])
  data
}

# ---------------------------------------------------------
# 10. QC reason prettifier
# ---------------------------------------------------------

prettify_reason_string <- function(x, dictionary = NULL) {
  if (is.null(dictionary)) {
    dictionary <- c(qc_flag_labels, exclusion_labels)
  }

  x <- as.character(x)

  prettified <- vapply(
    x,
    function(one_x) {
      if (is.na(one_x) || !nzchar(one_x)) {
        return(NA_character_)
      }

      pieces <- unlist(strsplit(one_x, ";\\s*"))
      pieces <- trimws(pieces)

      pieces <- ifelse(
        pieces %in% names(dictionary),
        dictionary[pieces],
        gsub("_", " ", pieces)
      )

      paste(unname(pieces), collapse = "; ")
    },
    character(1)
  )

  prettified
}

# ---------------------------------------------------------
# 11. Axis label convenience helpers
# ---------------------------------------------------------

xlab_kelp <- function(var_name, plain = FALSE) {
  ggplot2::xlab(label_variable(var_name, plain = plain))
}

ylab_kelp <- function(var_name, plain = FALSE) {
  ggplot2::ylab(label_variable(var_name, plain = plain))
}

labs_kelp <- function(title = NULL,
                      subtitle = NULL,
                      caption = NULL,
                      x = NULL,
                      y = NULL,
                      plain_axes = FALSE) {
  ggplot2::labs(
    title = title,
    subtitle = subtitle,
    caption = caption,
    x = if (!is.null(x)) label_variable(x, plain = plain_axes) else NULL,
    y = if (!is.null(y)) label_variable(y, plain = plain_axes) else NULL
  )
}

# ---------------------------------------------------------
# 12. Standard caption building helpers
# ---------------------------------------------------------

caption_model_ci <- function() {
  "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals."
}

caption_model_ci_day <- function() {
  "Points show raw observations; lines and ribbons show model predictions with 95% confidence intervals, facetted by day."
}

caption_qc_summary <- function() {
  "Summary of QC flags and explicit exclusions applied before analysis."
}

# ---------------------------------------------------------
# 13. Standard manuscript text snippets
# ---------------------------------------------------------

manuscript_snippets <- list(
  motility_definition = "Motile fraction was defined as mobile_cell_count divided by the sum of mobile_cell_count and stationary_cell_count.",
  primary_response = "Motility models were fitted from primary cell counts using cbind(mobile_cell_count, stationary_cell_count).",
  single_species_scope = "Manuscript 4 is restricted to Laminaria digitata; species is therefore not included as a model term."
)

get_snippet <- function(name) {
  if (!name %in% names(manuscript_snippets)) {
    stop(
      paste0(
        "Unknown snippet name: ", name,
        ". Available snippets are: ",
        paste(names(manuscript_snippets), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  manuscript_snippets[[name]]
}

# ---------------------------------------------------------
# 14. Validation helpers
# ---------------------------------------------------------

assert_known_labels <- function(x, dictionary_name = c(
  "experiment",
  "particle",
  "size_class",
  "season",
  "species",
  "variable",
  "qc_flag",
  "exclusion",
  "figure_title",
  "table_title"
)) {
  dictionary_name <- match.arg(dictionary_name)

  dictionary <- switch(
    dictionary_name,
    experiment   = experiment_labels,
    particle     = particle_labels,
    size_class   = size_class_labels,
    season       = season_labels,
    species      = species_labels,
    variable     = variable_labels_plain,
    qc_flag      = qc_flag_labels,
    exclusion    = exclusion_labels,
    figure_title = figure_title_labels,
    table_title  = table_title_labels
  )

  missing_values <- setdiff(unique(as.character(x)), names(dictionary))

  if (length(missing_values) > 0) {
    warning(
      paste0(
        "Unknown values for ", dictionary_name, " labels: ",
        paste(missing_values, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------
# 15. Message on load
# ---------------------------------------------------------

message("helpers_labels.R loaded successfully.")
message("Available dictionaries: experiment_labels, particle_labels, size_class_labels, variable_labels, qc_flag_labels")
message("Available helpers: label_experiment(), label_particle_type(), label_size_class(), label_variable(), prettify_reason_string(), labs_kelp()")

