# =========================================================
# Script title: 04_derive_variables_v3.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Derive analysis variables from the QC-flagged MS4
#          dataset for summaries, plotting, experiment-
#          specific modelling support, and downstream
#          figure/table generation.
# Inputs:
# - data_derived/ms4_qc_flagged_full.csv or .rds
# Outputs:
# - data_derived/ms4_analysis_derived.csv
# - data_derived/ms4_analysis_derived.rds
# - data_derived/ms4_analysis_derived.parquet (if arrow available)
# - outputs/tables/04_derived_variable_summary.csv
# - outputs/tables/04_experiment_level_summary.csv
# - outputs/tables/04_ms4_block_summary.csv
# - outputs/logs/04_derive_variables_log_*.txt
# Date created: 01 March 2026
# Last updated: 27 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 03_qc_and_exclusions.R before this script.
# - This script derives variables only; it does not apply
#   any new QC rules.
# - Primary modelling remains based on
#   cbind(mobile_cell_count, stationary_cell_count).
# - Manuscript 4 is restricted to the concurrent experiment
#   block 8.2, 9.2, 10.2, and 11.2.
# - experiment / experiment_num are the primary subsetting
#   keys; experiment_type is used for labelling/summaries.
# - For Experiment 3, two parallel modelling builds are
#   supported:
#     (a) main size-class comparison: BWC vs BWF only
#     (b) secondary brake-presence comparison:
#         CONTROL vs BRAKE_WEAR (BWC + BWF pooled)
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 04: DERIVE VARIABLES\n")
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
  "helpers_labels.R",
  "helpers_tables.R"
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

file_input_rds <- file.path(dir_data_derived, "ms4_qc_flagged_full.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_qc_flagged_full.csv")

file_output_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")
file_output_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_output_parquet <- file.path(dir_data_derived, "ms4_analysis_derived.parquet")

file_derived_summary <- file.path(dir_tables, "04_derived_variable_summary.csv")
file_experiment_summary <- file.path(dir_tables, "04_experiment_level_summary.csv")
file_ms4_block_summary <- file.path(dir_tables, "04_ms4_block_summary.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_derivation_log <- file.path(
  dir_logs,
  paste0("04_derive_variables_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 4. Import QC-flagged dataset
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
      "No QC-flagged dataset found.\nExpected one of:\n- ",
      file_input_rds,
      "\n- ",
      file_input_csv,
      "\nPlease run 03_qc_and_exclusions.R first."
    ),
    call. = FALSE
  )
}

cat("QC-flagged dataset loaded from:\n")
cat(input_source_used, "\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

# ---------------------------------------------------------
# 5. Check required columns
# ---------------------------------------------------------

required_columns <- c(
  "row_id",
  "experiment",
  "experiment_num",
  "experiment_type",
  "ms4_block_flag",
  "ms4_target_experiment",
  "mobile_cell_count",
  "stationary_cell_count",
  "days_from_start",
  "culture",
  "well",
  "exclusion_flag",
  "toxin_exposure"
)

missing_required_columns <- required_columns[!required_columns %in% names(dat)]

if (length(missing_required_columns) > 0) {
  stop(
    paste0(
      "The following required column(s) are missing:\n- ",
      paste(missing_required_columns, collapse = "\n- "),
      "\nPlease review upstream scripts."
    ),
    call. = FALSE
  )
}

cat("Required columns verified.\n\n")

# ---------------------------------------------------------
# 6. Restrict explicitly to the Manuscript 4 block
# ---------------------------------------------------------

ms4_experiment_nums <- c(8.2, 9.2, 10.2, 11.2)

dat <- dat |>
  dplyr::mutate(
    experiment_chr = stringr::str_squish(as.character(experiment)),
    experiment_num = suppressWarnings(as.numeric(as.character(experiment_num))),
    experiment_type = as.character(experiment_type),
    ms4_target_experiment = as.character(ms4_target_experiment),
    ms4_block_flag = as.logical(ms4_block_flag),
    exclusion_flag = as.logical(exclusion_flag)
  ) |>
  dplyr::filter(
    !is.na(experiment_num),
    experiment_num %in% ms4_experiment_nums
  )

cat("Restricted to Manuscript 4 experiment block.\n")
cat("Rows retained in MS4 block:", nrow(dat), "\n\n")

# ---------------------------------------------------------
# 7. Standardise key grouping fields
# ---------------------------------------------------------

standardise_upper_trim <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", ".")] <- NA_character_
  x
}

for (nm in intersect(
  c(
    "experiment_type",
    "ms4_target_experiment",
    "particle_type",
    "particle_class",
    "size_class",
    "toxin_exposure",
    "season",
    "site",
    "species",
    "species_std"
  ),
  names(dat)
)) {
  dat[[nm]] <- standardise_upper_trim(dat[[nm]])
}

cat("Key grouping fields standardised.\n\n")

# ---------------------------------------------------------
# 8. Refresh core count-derived variables
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    total_cells = dplyr::if_else(
      !is.na(mobile_cell_count) & !is.na(stationary_cell_count),
      as.integer(mobile_cell_count + stationary_cell_count),
      total_cells
    ),
    motility_ratio = dplyr::if_else(
      !is.na(mobile_cell_count) & !is.na(stationary_cell_count) &
        (mobile_cell_count + stationary_cell_count) > 0,
      mobile_cell_count / (mobile_cell_count + stationary_cell_count),
      NA_real_
    ),
    stationary_ratio = dplyr::if_else(
      !is.na(mobile_cell_count) & !is.na(stationary_cell_count) &
        (mobile_cell_count + stationary_cell_count) > 0,
      stationary_cell_count / (mobile_cell_count + stationary_cell_count),
      NA_real_
    )
  )

cat("Core count-derived variables refreshed.\n\n")

# ---------------------------------------------------------
# 9. Derive experiment identity helpers
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    experiment_num_chr = format(experiment_num, trim = TRUE, scientific = FALSE),
    experiment_label_core = dplyr::case_when(
      experiment_num == 9.2 ~ "LIGHT_ONLY",
      experiment_num == 10.2 ~ "DEFINED_PARTICLES",
      experiment_num == 11.2 ~ "BRAKE_SIZE",
      experiment_num == 8.2 ~ "FIELD_SPM",
      TRUE ~ experiment_type
    ),
    experiment_label_core = dplyr::coalesce(experiment_label_core, experiment_type),
    experiment_plot_group = dplyr::case_when(
      experiment_num == 9.2 ~ "EXP1_LIGHT",
      experiment_num == 10.2 ~ "EXP2_DEFINED_PARTICLES",
      experiment_num == 11.2 ~ "EXP3_BRAKE_SIZE",
      experiment_num == 8.2 ~ "EXP4_FIELD_SPM",
      TRUE ~ ms4_target_experiment
    )
  )

cat("Experiment identity helpers derived from experiment_num / experiment_type.\n\n")

# ---------------------------------------------------------
# 10. Derive transformed numeric predictors
# ---------------------------------------------------------

optional_numeric_cols <- intersect(
  c("lux_exposure", "ntu", "mass_loading_ug_l", "days_from_start"),
  names(dat)
)

for (nm in optional_numeric_cols) {
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
}

dat <- dat |>
  dplyr::mutate(
    log10_ntu_plus1 = dplyr::if_else(
      "ntu" %in% names(dat) & !is.na(ntu),
      log10(ntu + 1),
      NA_real_
    ),
    log10_mass_loading_ug_l_plus1 = dplyr::if_else(
      "mass_loading_ug_l" %in% names(dat) & !is.na(mass_loading_ug_l),
      log10(mass_loading_ug_l + 1),
      NA_real_
    ),
    lux_exposure_f = dplyr::case_when(
      !("lux_exposure" %in% names(dat)) ~ NA_character_,
      is.na(lux_exposure) ~ NA_character_,
      lux_exposure == 0 ~ "0 lux",
      lux_exposure == 4 ~ "4 lux",
      lux_exposure == 70 ~ "70 lux",
      lux_exposure == 117 ~ "117 lux",
      TRUE ~ paste0(lux_exposure, " lux")
    ),
    days_from_start_f = factor(
      days_from_start,
      levels = sort(unique(stats::na.omit(days_from_start)))
    ),
    day_order = dplyr::dense_rank(days_from_start)
  )

cat("Numeric transforms and day-factor variables derived.\n\n")

# ---------------------------------------------------------
# 11. Derive plotting and Experiment 3 comparison variables
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    particle_type_plot = dplyr::case_when(
      !is.na(particle_type) ~ particle_type,
      !is.na(toxin_exposure) ~ toxin_exposure,
      experiment_num == 8.2 ~ "SPM",
      TRUE ~ NA_character_
    ),
    particle_type_plot = dplyr::case_when(
      particle_type_plot %in% c("BRAKE WEAR COARSE", "BRAKE_WEAR_COARSE", "COARSE", "BWC") ~ "BWC",
      particle_type_plot %in% c("BRAKE WEAR FINE", "BRAKE_WEAR_FINE", "FINE", "BWF") ~ "BWF",
      particle_type_plot %in% c("BRAKE_WEAR", "BRAKE WEAR") ~ "BRAKE_WEAR",
      particle_type_plot %in% c("FIELD SPM", "FIELD_SPM", "SPM") ~ "SPM",
      TRUE ~ particle_type_plot
    ),
    size_class_plot = dplyr::case_when(
      !is.na(size_class) & size_class %in% c("COARSE", "BWC") ~ "COARSE",
      !is.na(size_class) & size_class %in% c("FINE", "BWF") ~ "FINE",
      particle_type_plot == "BWC" ~ "COARSE",
      particle_type_plot == "BWF" ~ "FINE",
      TRUE ~ NA_character_
    ),
    treatment_group = dplyr::case_when(
      experiment_num == 9.2 ~ "LIGHT_ONLY",
      experiment_num == 10.2 ~ "DEFINED_PARTICLES",
      experiment_num == 11.2 ~ "BRAKE_SIZE",
      experiment_num == 8.2 ~ "FIELD_SPM",
      TRUE ~ "OTHER"
    ),
    brake_exposure = dplyr::case_when(
      experiment_num == 11.2 & toxin_exposure %in% c("BWC", "BWF") ~ "BRAKE_WEAR",
      experiment_num == 11.2 & toxin_exposure == "CONTROL" ~ "CONTROL",
      TRUE ~ NA_character_
    ),
    brake_size_treatment = dplyr::case_when(
      experiment_num == 11.2 & toxin_exposure == "BWC" ~ "BWC",
      experiment_num == 11.2 & toxin_exposure == "BWF" ~ "BWF",
      TRUE ~ NA_character_
    )
  )

cat("Treatment, plotting, and Experiment 3 comparison variables derived.\n\n")

# ---------------------------------------------------------
# 12. Derive experiment-specific inclusion markers
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    use_exp1_light =
      experiment_num == 9.2 &
      !exclusion_flag &
      !is.na(lux_exposure) &
      !is.na(days_from_start),

    use_exp2_defined_particles =
      experiment_num == 10.2 &
      !exclusion_flag &
      !is.na(ntu) &
      !is.na(particle_type_plot) &
      !is.na(days_from_start),

    use_exp3_brake_size =
      experiment_num == 11.2 &
      !exclusion_flag &
      toxin_exposure %in% c("BWC", "BWF") &
      !is.na(size_class_plot) &
      !is.na(days_from_start),

    use_exp3_brake_presence =
      experiment_num == 11.2 &
      toxin_exposure %in% c("BWC", "BWF", "CONTROL") &
      !is.na(days_from_start) &
      !is.na(brake_exposure),

    use_exp4_field_spm =
      experiment_num == 8.2 &
      !exclusion_flag &
      !is.na(ntu) &
      !is.na(days_from_start)
  )

cat("Experiment-specific inclusion markers derived.\n\n")

# ---------------------------------------------------------
# 13. Derive plotting / faceting helper variables
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    day_label = dplyr::case_when(
      is.na(days_from_start) ~ NA_character_,
      TRUE ~ paste0("Day ", days_from_start)
    ),
    ntu_band = dplyr::case_when(
      !("ntu" %in% names(dat)) ~ NA_character_,
      is.na(ntu) ~ NA_character_,
      ntu < 25 ~ "<25 NTU",
      ntu >= 25 & ntu < 100 ~ "25-<100 NTU",
      ntu >= 100 & ntu < 400 ~ "100-<400 NTU",
      ntu >= 400 ~ ">=400 NTU",
      TRUE ~ NA_character_
    ),
    lux_band = dplyr::case_when(
      !("lux_exposure" %in% names(dat)) ~ NA_character_,
      is.na(lux_exposure) ~ NA_character_,
      lux_exposure == 0 ~ "0 lux",
      lux_exposure == 4 ~ "4 lux",
      lux_exposure == 70 ~ "70 lux",
      lux_exposure == 117 ~ "117 lux",
      TRUE ~ paste0(lux_exposure, " lux")
    )
  )

cat("Plotting and faceting helpers derived.\n\n")

# ---------------------------------------------------------
# 14. Derive within-experiment observation counts
# ---------------------------------------------------------

dat <- dat |>
  dplyr::group_by(experiment_num, days_from_start) |>
  dplyr::mutate(
    n_obs_experiment_day = dplyr::n()
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(experiment_num, culture) |>
  dplyr::mutate(
    n_obs_experiment_culture = dplyr::n()
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(experiment_num, well) |>
  dplyr::mutate(
    n_obs_experiment_well = dplyr::n()
  ) |>
  dplyr::ungroup()

cat("Within-experiment observation counts derived.\n\n")

# ---------------------------------------------------------
# 15. Derive human-readable labels
# ---------------------------------------------------------

if (exists("label_experiment", mode = "function", inherits = TRUE)) {
  dat$experiment_label <- label_experiment(dat$experiment_label_core, short = FALSE)
  dat$experiment_short_label <- label_experiment(dat$experiment_label_core, short = TRUE)
} else {
  dat$experiment_label <- dat$experiment_plot_group
  dat$experiment_short_label <- dat$experiment_plot_group
}

if (exists("label_particle_type", mode = "function", inherits = TRUE)) {
  dat$particle_label <- label_particle_type(dat$particle_type_plot, short = FALSE)
  dat$particle_short_label <- label_particle_type(dat$particle_type_plot, short = TRUE)
} else {
  dat$particle_label <- dat$particle_type_plot
  dat$particle_short_label <- dat$particle_type_plot
}

if (exists("label_size_class", mode = "function", inherits = TRUE)) {
  dat$size_class_label <- label_size_class(dat$size_class_plot)
} else {
  dat$size_class_label <- dat$size_class_plot
}

cat("Human-readable labels derived.\n\n")

# ---------------------------------------------------------
# 16. Re-coerce key grouping variables as factors
# ---------------------------------------------------------

factor_cols <- intersect(
  c(
    "experiment",
    "experiment_num",
    "experiment_num_chr",
    "experiment_type",
    "experiment_label_core",
    "experiment_plot_group",
    "experiment_label",
    "experiment_short_label",
    "ms4_target_experiment",
    "treatment_group",
    "particle_type_plot",
    "particle_label",
    "particle_short_label",
    "size_class_plot",
    "size_class_label",
    "brake_exposure",
    "brake_size_treatment",
    "day_label",
    "days_from_start_f",
    "ntu_band",
    "lux_band",
    "lux_exposure_f",
    "culture",
    "well",
    "species",
    "species_std",
    "site",
    "season",
    "collection_site",
    "toxin_exposure"
  ),
  names(dat)
)

for (nm in factor_cols) {
  dat[[nm]] <- as.factor(dat[[nm]])
}

cat("Key grouping variables re-coerced as factors.\n\n")

# ---------------------------------------------------------
# 17. Build derived variable summary table
# ---------------------------------------------------------

derived_variable_summary <- tibble::tibble(
  derived_variable = c(
    "total_cells",
    "motility_ratio",
    "stationary_ratio",
    "experiment_label_core",
    "experiment_plot_group",
    "log10_ntu_plus1",
    "log10_mass_loading_ug_l_plus1",
    "lux_exposure_f",
    "days_from_start_f",
    "day_order",
    "treatment_group",
    "particle_type_plot",
    "size_class_plot",
    "brake_exposure",
    "brake_size_treatment",
    "use_exp1_light",
    "use_exp2_defined_particles",
    "use_exp3_brake_size",
    "use_exp3_brake_presence",
    "use_exp4_field_spm",
    "day_label",
    "ntu_band",
    "lux_band",
    "n_obs_experiment_day",
    "n_obs_experiment_culture",
    "n_obs_experiment_well",
    "experiment_label",
    "experiment_short_label",
    "particle_label",
    "particle_short_label",
    "size_class_label"
  ),
  description = c(
    "Recalculated mobile + stationary count total",
    "Mobile / total; for summaries and plotting only",
    "Stationary / total; complementary descriptive ratio",
    "Experiment label derived from experiment_num",
    "MS4 experiment plotting label",
    "log10(ntu + 1) transform for NTU-based analyses",
    "log10(mass_loading_ug_l + 1) transform for mass-based analyses",
    "Factor version of lux_exposure for plotting/sensitivity work",
    "Factor version of days_from_start",
    "Ordered integer rank of day values",
    "Top-level experiment grouping variable",
    "Plotting treatment factor for particle exposures",
    "Plotting size factor for brake-size comparison",
    "Secondary Experiment 3 predictor: any brake wear vs control",
    "Primary Experiment 3 treatment identity: BWC vs BWF",
    "Marker for rows eligible for Experiment 1 analysis",
    "Marker for rows eligible for Experiment 2 analysis",
    "Marker for rows eligible for main Experiment 3 size-class analysis",
    "Marker for rows eligible for secondary Experiment 3 brake-presence analysis",
    "Marker for rows eligible for Experiment 4 analysis",
    "Human-readable day label for faceting",
    "Broad NTU bin for summaries/EDA",
    "Human-readable irradiance level label",
    "Observation count within experiment-day",
    "Observation count within experiment-culture",
    "Observation count within experiment-well",
    "Full experiment label",
    "Short experiment label",
    "Full particle/treatment label",
    "Short particle/treatment label",
    "Human-readable size class label"
  )
)

readr::write_csv(derived_variable_summary, file_derived_summary)

cat("Derived variable summary written to:\n")
cat(file_derived_summary, "\n\n")

# ---------------------------------------------------------
# 18. Build experiment-level summary table
# ---------------------------------------------------------

safe_min <- function(x) {
  if (!is.numeric(x)) return(NA_real_)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  if (!is.numeric(x)) return(NA_real_)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

experiment_level_summary <- dat |>
  dplyr::group_by(
    experiment_num,
    experiment_label,
    days_from_start
  ) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_rows_not_excluded = sum(!exclusion_flag, na.rm = TRUE),
    n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
    n_wells = dplyr::n_distinct(well, na.rm = TRUE),
    n_particles = dplyr::n_distinct(particle_type_plot, na.rm = TRUE),
    mean_total_cells = mean(total_cells, na.rm = TRUE),
    mean_motility_ratio = mean(motility_ratio, na.rm = TRUE),
    min_ntu = safe_min(ntu),
    max_ntu = safe_max(ntu),
    min_lux = safe_min(lux_exposure),
    max_lux = safe_max(lux_exposure),
    min_mass_loading_ug_l = safe_min(mass_loading_ug_l),
    max_mass_loading_ug_l = safe_max(mass_loading_ug_l),
    .groups = "drop"
  ) |>
  dplyr::arrange(experiment_num, days_from_start)

readr::write_csv(experiment_level_summary, file_experiment_summary)

cat("Experiment-level summary written to:\n")
cat(file_experiment_summary, "\n\n")

# ---------------------------------------------------------
# 19. Build MS4 block summary table
# ---------------------------------------------------------

has_video_file <- "video_file" %in% names(dat)

ms4_block_summary <- dat |>
  dplyr::group_by(
    experiment_num,
    experiment_label,
    experiment_short_label
  ) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_rows_not_excluded = sum(!exclusion_flag, na.rm = TRUE),
    n_video_files = if (has_video_file) dplyr::n_distinct(video_file, na.rm = TRUE) else NA_integer_,
    n_cultures = dplyr::n_distinct(culture, na.rm = TRUE),
    n_wells = dplyr::n_distinct(well, na.rm = TRUE),
    n_days = dplyr::n_distinct(days_from_start, na.rm = TRUE),
    n_use_exp1_light = sum(use_exp1_light, na.rm = TRUE),
    n_use_exp2_defined_particles = sum(use_exp2_defined_particles, na.rm = TRUE),
    n_use_exp3_brake_size = sum(use_exp3_brake_size, na.rm = TRUE),
    n_use_exp3_brake_presence = sum(use_exp3_brake_presence, na.rm = TRUE),
    n_use_exp4_field_spm = sum(use_exp4_field_spm, na.rm = TRUE),
    min_day = safe_min(days_from_start),
    max_day = safe_max(days_from_start),
    mean_total_cells = mean(total_cells, na.rm = TRUE),
    mean_motility_ratio = mean(motility_ratio, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(experiment_num)

readr::write_csv(ms4_block_summary, file_ms4_block_summary)

cat("MS4 block summary written to:\n")
cat(file_ms4_block_summary, "\n\n")

# ---------------------------------------------------------
# 20. Write derived dataset to disk
# ---------------------------------------------------------

readr::write_csv(dat, file_output_csv)
saveRDS(dat, file_output_rds)

cat("Derived dataset written to:\n")
cat("- ", file_output_csv, "\n", sep = "")
cat("- ", file_output_rds, "\n", sep = "")

if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(dat, file_output_parquet)
  cat("- ", file_output_parquet, "\n", sep = "")
} else {
  cat("Parquet output skipped because package 'arrow' is not available.\n")
}

cat("\n")

# ---------------------------------------------------------
# 21. Write derivation log
# ---------------------------------------------------------

sink(file_derivation_log)
cat("SPM Analysis - 04_derive_variables log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")

cat("Input dataset:\n")
cat(input_source_used, "\n\n")

cat("Output files:\n")
cat(file_output_csv, "\n")
cat(file_output_rds, "\n")
if (requireNamespace("arrow", quietly = TRUE)) {
  cat(file_output_parquet, "\n")
}
cat(file_derived_summary, "\n")
cat(file_experiment_summary, "\n")
cat(file_ms4_block_summary, "\n\n")

cat("Dimensions:\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

cat("Derived variable summary:\n")
print(derived_variable_summary)

cat("\nExperiment-level summary:\n")
print(experiment_level_summary)

cat("\nMS4 block summary:\n")
print(ms4_block_summary)

cat("\nSession information:\n\n")
print(utils::sessionInfo())
sink()

cat("Derivation log written to:\n")
cat(file_derivation_log, "\n\n")

# ---------------------------------------------------------
# 22. Console summary
# ---------------------------------------------------------

cat("Experiment eligibility counts:\n")
print(
  c(
    use_exp1_light = sum(dat$use_exp1_light, na.rm = TRUE),
    use_exp2_defined_particles = sum(dat$use_exp2_defined_particles, na.rm = TRUE),
    use_exp3_brake_size = sum(dat$use_exp3_brake_size, na.rm = TRUE),
    use_exp3_brake_presence = sum(dat$use_exp3_brake_presence, na.rm = TRUE),
    use_exp4_field_spm = sum(dat$use_exp4_field_spm, na.rm = TRUE)
  )
)
cat("\n")

cat("Experiment counts by experiment_num:\n")
print(table(dat$experiment_num, useNA = "ifany"))
cat("\n")

cat("Experiment counts by experiment_type:\n")
print(table(dat$experiment_type, useNA = "ifany"))
cat("\n")

cat("Experiment 3 size-class rows by toxin_exposure:\n")
print(table(dat$toxin_exposure[dat$use_exp3_brake_size], useNA = "ifany"))
cat("\n")

cat("Experiment 3 brake-presence rows by toxin_exposure:\n")
print(table(dat$toxin_exposure[dat$use_exp3_brake_presence], useNA = "ifany"))
cat("\n")

cat("========================================================\n")
cat("SCRIPT 04 COMPLETE: DERIVE VARIABLES\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

