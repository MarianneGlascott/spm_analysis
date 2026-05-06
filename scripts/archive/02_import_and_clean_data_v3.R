# =========================================================
# Script title: 02_import_and_clean_data_v5.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Import the raw master dataset, standardise column
#          names and data types, retain the variables needed
#          for Manuscript 4, derive key manuscript flags and
#          experiment labels, and write a clean imported
#          dataset to disk without applying exclusions.
# Inputs: data_raw/All_data_cleaned.csv
# Outputs:
# - data_clean/ms4_imported_clean.csv
# - data_clean/ms4_imported_clean.rds
# - data_clean/ms4_imported_clean.parquet (if arrow available)
# - outputs/logs/02_import_and_clean_data_log_*.txt
# - outputs/tables/02_variable_classes_import.csv
# - outputs/tables/02_experiment_code_summary.csv
# Date created: 24 March 2026
# Last updated: 26 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - This script does not apply exclusions.
# - Raw data must remain untouched.
# - The master dataset contains experiments beyond those used
#   in Manuscript 4.
# - Manuscript 4 experiment block:
#   8.2 = field SPM
#   9.2 = light-only
#   10.2 = defined particles
#   11.2 = brake-wear size comparison
# - QC flags and exclusions are handled in
#   03_qc_and_exclusions.R
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 02: IMPORT AND CLEAN DATA\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

source(here::here("R", "helpers_theme.R"))
source(here::here("R", "helpers_save_figures.R"))
source(here::here("R", "helpers_labels.R"))
source(here::here("R", "helpers_tables.R"))
source(here::here("R", "helpers_model_checks.R"))

# ---------------------------------------------------------
# 1. Check that setup objects exist
# ---------------------------------------------------------

required_objects <- c(
  "project_root",
  "dir_data_clean",
  "dir_logs",
  "dir_tables",
  "file_raw_master_csv",
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
# 2. Check source file exists
# ---------------------------------------------------------

if (!file.exists(file_raw_master_csv)) {
  stop(
    paste0(
      "Raw data file not found at expected location:\n",
      file_raw_master_csv,
      "\nPlease check the file path and filename."
    ),
    call. = FALSE
  )
}

cat("Raw source file found:\n")
cat(file_raw_master_csv, "\n\n")

# ---------------------------------------------------------
# 3. Define output paths
# ---------------------------------------------------------

file_clean_csv <- file.path(dir_data_clean, "ms4_imported_clean.csv")
file_clean_rds <- file.path(dir_data_clean, "ms4_imported_clean.rds")
file_clean_parquet <- file.path(dir_data_clean, "ms4_imported_clean.parquet")

file_variable_classes <- file.path(dir_tables, "02_variable_classes_import.csv")
file_experiment_code_summary <- file.path(dir_tables, "02_experiment_code_summary.csv")

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_import_log <- file.path(
  dir_logs,
  paste0("02_import_and_clean_data_log_", timestamp_now, ".txt")
)

# ---------------------------------------------------------
# 4. Import raw data
# ---------------------------------------------------------

cat("Importing raw CSV...\n")

raw_dat <- readr::read_csv(
  file = file_raw_master_csv,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "N/A", "NULL", "null", ".", " ")
)

cat("Raw data imported successfully.\n")
cat("Rows:", nrow(raw_dat), "\n")
cat("Columns:", ncol(raw_dat), "\n\n")

# ---------------------------------------------------------
# 5. Preserve original names and standardise column names
# ---------------------------------------------------------

original_names <- names(raw_dat)

dat <- raw_dat |>
  janitor::clean_names()

cat("Column names standardised with janitor::clean_names().\n\n")

name_key <- tibble::tibble(
  original_name = original_names,
  cleaned_name = names(dat)
)

# ---------------------------------------------------------
# 6. Add stable row identifier
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(row_id = dplyr::row_number(), .before = 1)

cat("Stable row_id created.\n\n")

# ---------------------------------------------------------
# 7. Harmonise likely alternative column names
# ---------------------------------------------------------

rename_if_present <- function(data, old, new) {
  if (old %in% names(data) && !new %in% names(data)) {
    data <- dplyr::rename(data, !!new := !!rlang::sym(old))
  }
  data
}

dat <- dat |>
  rename_if_present("experiment_type_name", "experiment_type") |>
  rename_if_present("count_date", "count_date_yyyy_mm_dd") |>
  rename_if_present("start_date", "start_date_yyyy_mm_dd") |>
  rename_if_present("video_file_name", "video_file") |>
  rename_if_present("cama_file_name", "camera_file")

cat("Common alternative names harmonised where present.\n\n")

# ---------------------------------------------------------
# 8. Ensure expected core columns exist where possible
# ---------------------------------------------------------

expected_core_columns <- c(
  "row_id",
  "video_file",
  "camera_file",
  "frame_count",
  "mobile_cell_count",
  "stationary_cell_count",
  "total_cells",
  "motility_ratio",
  "count_date_yyyy_mm_dd",
  "start_date_yyyy_mm_dd",
  "species",
  "site",
  "culture",
  "experiment",
  "experiment_id",
  "experiment_num",
  "experiment_type",
  "well",
  "days_from_start",
  "toxin_exposure",
  "particle_type",
  "particle_class",
  "size_class",
  "cu_ug_l",
  "bw_ug_l",
  "ntu",
  "tile",
  "tile_count",
  "lux_exposure",
  "mass_loading_ug_l",
  "comment",
  "collection_date",
  "year",
  "season",
  "collection_site",
  "preparation",
  "release"
)

missing_expected_columns <- setdiff(expected_core_columns, names(dat))

if (length(missing_expected_columns) > 0) {
  cat("Expected columns not currently present in source data:\n")
  cat("-", paste(missing_expected_columns, collapse = "\n- "), "\n\n")
  
  for (nm in missing_expected_columns) {
    dat[[nm]] <- NA
  }
  
  cat("Missing expected columns added as NA placeholders.\n\n")
} else {
  cat("All expected core columns are present.\n\n")
}
# ---------------------------------------------------------
# 9. Coerce dates where present
# ---------------------------------------------------------

date_cols <- intersect(
  c(
    "count_date_yyyy_mm_dd",
    "start_date_yyyy_mm_dd",
    "collection_date",
    "preparation",
    "release"
  ),
  names(dat)
)

for (nm in date_cols) {
  dat[[nm]] <- suppressWarnings(lubridate::ymd(dat[[nm]]))
}

cat("Date columns coerced where present.\n\n")

# ---------------------------------------------------------
# 10. Coerce numeric/integer columns where present
# ---------------------------------------------------------

integer_cols <- intersect(
  c(
    "frame_count",
    "mobile_cell_count",
    "stationary_cell_count",
    "total_cells",
    "tile_count",
    "days_from_start"
  ),
  names(dat)
)

numeric_cols <- intersect(
  c(
    "motility_ratio",
    "cu_ug_l",
    "bw_ug_l",
    "ntu",
    "lux_exposure",
    "mass_loading_ug_l"
  ),
  names(dat)
)

for (nm in integer_cols) {
  dat[[nm]] <- suppressWarnings(as.integer(dat[[nm]]))
}

for (nm in numeric_cols) {
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
}

cat("Numeric and integer columns coerced where present.\n\n")

# ---------------------------------------------------------
# 11. Coerce grouping / metadata variables
# ---------------------------------------------------------

factor_cols <- intersect(
  c(
    "video_file",
    "camera_file",
    "species",
    "site",
    "culture",
    "experiment",
    "experiment_id",
    "experiment_type",
    "well",
    "toxin_exposure",
    "particle_type",
    "particle_class",
    "size_class",
    "tile",
    "year",
    "season",
    "collection_site"
  ),
  names(dat)
)

for (nm in factor_cols) {
  dat[[nm]] <- as.factor(dat[[nm]])
}

cat("Grouping and metadata columns coerced to factor where present.\n\n")

# ---------------------------------------------------------
# 12. Refresh total_cells and motility_ratio from primary counts
# ---------------------------------------------------------

if (all(c("mobile_cell_count", "stationary_cell_count") %in% names(dat))) {
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
        motility_ratio
      )
    )
}

cat("total_cells and motility_ratio refreshed where primary counts were available.\n\n")

# ---------------------------------------------------------
# 13. Derive experiment code helpers and manuscript flags
# ---------------------------------------------------------

dat <- dat |>
  dplyr::mutate(
    experiment_chr = as.character(experiment),
    experiment_chr = stringr::str_squish(experiment_chr),
    experiment_chr = dplyr::na_if(experiment_chr, ""),
    experiment_num = suppressWarnings(as.numeric(experiment_chr)),
    
    experiment_type = dplyr::case_when(
      experiment_num == 9.2  ~ "LIGHT_ONLY",
      experiment_num == 10.2 ~ "DEFINED_PARTICLES",
      experiment_num == 11.2 ~ "BRAKE_SIZE",
      experiment_num == 8.2  ~ "FIELD_SPM",
      TRUE ~ as.character(experiment_type)
    ),
    experiment_type = dplyr::na_if(experiment_type, ""),
    
    ms4_block_flag = !is.na(experiment_num) & experiment_num %in% c(8.2, 9.2, 10.2, 11.2),
    
    ms4_target_experiment = dplyr::case_when(
      experiment_num == 9.2  ~ "EXP1_LIGHT",
      experiment_num == 10.2 ~ "EXP2_DEFINED_PARTICLES",
      experiment_num == 11.2 ~ "EXP3_BRAKE_SIZE",
      experiment_num == 8.2  ~ "EXP4_FIELD_SPM",
      TRUE ~ NA_character_
    ),
    
    species_std = dplyr::case_when(
      is.na(species) ~ NA_character_,
      stringr::str_squish(as.character(species)) == "L. digitata" ~ "Laminaria digitata",
      stringr::str_squish(as.character(species)) == "L. hyperborea" ~ "Laminaria hyperborea",
      stringr::str_squish(as.character(species)) == "S. latissima" ~ "Saccharina latissima",
      TRUE ~ stringr::str_squish(as.character(species))
    ),
    
    particle_type = dplyr::case_when(
      experiment_type == "DEFINED_PARTICLES" & as.character(toxin_exposure) %in% c("SAND", "KAOLINITE", "PEAT") ~ as.character(toxin_exposure),
      experiment_type == "FIELD_SPM" & as.character(toxin_exposure) == "SPM" ~ "SPM",
      experiment_type == "BRAKE_SIZE" & as.character(toxin_exposure) %in% c("BWC", "BWF") ~ "BRAKE_WEAR",
      TRUE ~ as.character(particle_type)
    ),
    particle_type = dplyr::na_if(stringr::str_squish(particle_type), ""),
    
    size_class = dplyr::case_when(
      experiment_type == "BRAKE_SIZE" & as.character(toxin_exposure) == "BWC" ~ "COARSE",
      experiment_type == "BRAKE_SIZE" & as.character(toxin_exposure) == "BWF" ~ "FINE",
      TRUE ~ as.character(size_class)
    ),
    size_class = dplyr::na_if(stringr::str_squish(size_class), ""),
    
    is_control = !is.na(toxin_exposure) & as.character(toxin_exposure) == "CONTROL",
    
    is_target_species = dplyr::case_when(
      is.na(species_std) ~ NA,
      species_std == "Laminaria digitata" ~ TRUE,
      TRUE ~ FALSE
    )
  )

dat$experiment_type <- as.factor(dat$experiment_type)
dat$ms4_target_experiment <- as.factor(dat$ms4_target_experiment)
dat$species_std <- as.factor(dat$species_std)
dat$particle_type <- as.factor(dat$particle_type)
dat$size_class <- as.factor(dat$size_class)

cat("Experiment helpers and manuscript flags derived.\n")
cat("Note: no experiment filtering has been applied at import stage.\n\n")

# ---------------------------------------------------------
# 13b. Derive brake-wear mass loading (bw_ug_l)
# ---------------------------------------------------------

if ("cu_ug_l" %in% names(dat)) {
  
  dat <- dat |>
    dplyr::mutate(
      bw_ug_l = dplyr::case_when(
        experiment_num == 11.2 ~ cu_ug_l,
        TRUE ~ NA_real_
      )
    )
  
  cat("bw_ug_l derived from cu_ug_l for Experiment 3.\n\n")
}
# ---------------------------------------------------------
# 14. Build experiment code summary table
# ---------------------------------------------------------

experiment_code_summary <- dat |>
  dplyr::mutate(
    experiment_display = dplyr::case_when(
      is.na(experiment_chr) ~ "Missing",
      TRUE ~ experiment_chr
    )
  ) |>
  dplyr::group_by(
    experiment_display,
    experiment_num,
    experiment_type,
    ms4_block_flag,
    ms4_target_experiment
  ) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_controls = sum(is_control, na.rm = TRUE),
    n_target_species = sum(is_target_species, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(ms4_block_flag), experiment_num)

readr::write_csv(experiment_code_summary, file_experiment_code_summary)

cat("Experiment code summary written to:\n")
cat(file_experiment_code_summary, "\n\n")

# ---------------------------------------------------------
# 15. Keep all available columns, but order key analysis vars
# ---------------------------------------------------------

priority_cols <- c(
  "row_id",
  "video_file",
  "camera_file",
  "frame_count",
  "experiment",
  "experiment_id",
  "experiment_chr",
  "experiment_num",
  "ms4_block_flag",
  "ms4_target_experiment",
  "experiment_type",
  "is_control",
  "is_target_species",
  "culture",
  "well",
  "days_from_start",
  "species",
  "site",
  "season",
  "collection_site",
  "toxin_exposure",
  "particle_type",
  "particle_class",
  "size_class",
  "lux_exposure",
  "ntu",
  "mass_loading_ug_l",
  "cu_ug_l",
  "bw_ug_l",
  "mobile_cell_count",
  "stationary_cell_count",
  "total_cells",
  "motility_ratio",
  "tile",
  "tile_count",
  "count_date_yyyy_mm_dd",
  "start_date_yyyy_mm_dd",
  "collection_date",
  "release",
  "year",
  "preparation",
  "comment"
)

priority_cols <- intersect(priority_cols, names(dat))
other_cols <- setdiff(names(dat), priority_cols)

dat <- dat[, c(priority_cols, other_cols)]

cat("Column order finalised.\n\n")

# ---------------------------------------------------------
# 16. Create variable classes summary
# ---------------------------------------------------------

variable_classes <- tibble::tibble(
  column_name = names(dat),
  class = vapply(dat, function(x) class(x)[1], character(1)),
  n_missing = vapply(dat, function(x) sum(is.na(x)), integer(1)),
  n_unique = vapply(dat, function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1))
)

readr::write_csv(variable_classes, file_variable_classes)

cat("Variable class summary written to:\n")
cat(file_variable_classes, "\n\n")

# ---------------------------------------------------------
# 17. Write cleaned imported data to disk
# ---------------------------------------------------------

readr::write_csv(dat, file_clean_csv)
saveRDS(dat, file_clean_rds)

cat("Clean imported dataset written to:\n")
cat("- ", file_clean_csv, "\n", sep = "")
cat("- ", file_clean_rds, "\n", sep = "")

if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(dat, file_clean_parquet)
  cat("- ", file_clean_parquet, "\n", sep = "")
} else {
  cat("Parquet output skipped because package 'arrow' is not available.\n")
}

cat("\n")

# ---------------------------------------------------------
# 18. Write import log
# ---------------------------------------------------------

sink(file_import_log)
cat("SPM Analysis - 02_import_and_clean_data log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")

cat("Source file:\n")
cat(file_raw_master_csv, "\n\n")

cat("Output files:\n")
cat(file_clean_csv, "\n")
cat(file_clean_rds, "\n")
if (requireNamespace("arrow", quietly = TRUE)) {
  cat(file_clean_parquet, "\n")
}
cat("\n")

cat("Variable classes summary file:\n")
cat(file_variable_classes, "\n\n")

cat("Experiment code summary file:\n")
cat(file_experiment_code_summary, "\n\n")

cat("Dataset dimensions:\n")
cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

cat("Missing expected columns at import stage:\n")
if (length(missing_expected_columns) > 0) {
  cat(paste(missing_expected_columns, collapse = "\n"), "\n\n")
} else {
  cat("None\n\n")
}

cat("Experiment code summary:\n")
print(experiment_code_summary)
cat("\n")

cat("Variable classes:\n")
print(variable_classes)
cat("\n")

cat("Session information:\n\n")
print(utils::sessionInfo())
sink()

cat("Import log written to:\n")
cat(file_import_log, "\n\n")

# ---------------------------------------------------------
# 19. Console summaries
# ---------------------------------------------------------

cat("Top experiment summary rows:\n")
print(utils::head(experiment_code_summary, 10))
cat("\n")

cat("Top variable class summary rows:\n")
print(utils::head(variable_classes, 10))
cat("\n")

cat("========================================================\n")
cat("SCRIPT 02 COMPLETE: IMPORT AND CLEAN DATA\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

