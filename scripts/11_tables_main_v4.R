# =========================================================
# Script title: 11_tables_main_v4.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Generate the main publication tables for
#          Manuscript 4 from the derived dataset and saved
#          experiment-specific model outputs.
# Inputs:
# - data_derived/ms4_analysis_derived.csv or .rds
# - outputs/tables/06_exp1_light_dataset_summary.csv
# - outputs/tables/07_exp2_defined_particles_dataset_summary.csv
# - outputs/tables/08_exp3_brake_size_dataset_summary.csv
# - outputs/tables/09_exp4_field_spm_dataset_summary.csv
# - outputs/tables/06_exp1_light_model_summary_table.csv
# - outputs/tables/07_exp2_defined_particles_model_summary_table.csv
# - outputs/tables/08_exp3_brake_size_model_summary_table.csv
# - outputs/tables/09_exp4_field_spm_model_summary_table.csv
# - outputs/tables/06_exp1_light_model_comparison.csv
# - outputs/tables/07_exp2_defined_particles_model_comparison.csv
# - outputs/tables/08_exp3_brake_size_model_comparison.csv
# - outputs/tables/09_exp4_field_spm_model_comparison.csv
# - outputs/models/exp4_field_spm/exp4_field_spm_day4_model_comparison.csv
# Outputs:
# - outputs/tables/Table1_experimental_design_summary.{csv,pdf,png,tiff}
# - outputs/tables/Table2_dataset_summary.{csv,pdf,png,tiff}
# - outputs/tables/Table3_model_summary_table.{csv,pdf,png,tiff}
# - outputs/tables/Table4_model_comparison_summary.{csv,pdf,png,tiff}
# - outputs/tables/11_tables_main_manifest.csv
# - outputs/logs/11_tables_main_log_*.txt
# Date created: 26 March 2026
# Last updated: 31 March 2026
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 11: TABLES MAIN\n")
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
# 2. Package checks
# ---------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "tibble", "ggplot2", "gridExtra")

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
# 3. Small local helpers
# ---------------------------------------------------------

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

write_caption_md_local <- function(table_name, caption_text) {
  dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
  caption_file <- file.path(dir_tables, paste0(table_name, "_caption.md"))
  
  writeLines(
    c(
      paste0("# ", table_name),
      "",
      caption_text
    ),
    con = caption_file
  )
  
  invisible(caption_file)
}

save_table_multi <- function(data,
                             file_stem,
                             caption_text = NULL,
                             width = 11,
                             height = NULL,
                             base_size = 10) {
  dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
  
  file_csv <- file.path(dir_tables, paste0(file_stem, ".csv"))
  file_pdf <- file.path(dir_tables, paste0(file_stem, ".pdf"))
  file_png <- file.path(dir_tables, paste0(file_stem, ".png"))
  file_tiff <- file.path(dir_tables, paste0(file_stem, ".tiff"))
  
  readr::write_csv(data, file_csv)
  
  n_rows <- nrow(data)
  if (is.null(height)) {
    height <- max(2.5, 1.2 + 0.35 * (n_rows + 1))
  }
  
  tg <- gridExtra::tableGrob(
    data,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = base_size,
      core = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  
  grDevices::cairo_pdf(file_pdf, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(tg)
  grDevices::dev.off()
  
  ggplot2::ggsave(
    filename = file_png,
    plot = tg,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file_tiff,
    plot = tg,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  caption_file <- NULL
  if (!is.null(caption_text) && nzchar(caption_text)) {
    caption_file <- write_caption_md_local(file_stem, caption_text)
  }
  
  tibble::tibble(
    table_name = file_stem,
    file_csv = file_csv,
    file_pdf = file_pdf,
    file_png = file_png,
    file_tiff = file_tiff,
    caption_file = caption_file
  )
}

extract_h4 <- function(df, label) {
  if (is.null(df)) return(NULL)
  
  mono <- df |>
    dplyr::filter(model == "monotonic_primary") |>
    dplyr::pull(AIC)
  
  nonlin <- df |>
    dplyr::filter(model == "nonlinear_primary") |>
    dplyr::pull(AIC)
  
  if (length(mono) == 0 || length(nonlin) == 0) return(NULL)
  
  delta_aic <- mono[1] - nonlin[1]
  
  tibble::tibble(
    dataset = label,
    monotonic_aic = mono[1],
    nonlinear_aic = nonlin[1],
    delta_aic = delta_aic,
    interpretation = dplyr::case_when(
      delta_aic > 10 ~ "Strong support for nonlinear NTU-response",
      delta_aic > 2 ~ "Moderate support for nonlinear NTU-response",
      TRUE ~ "No strong support for nonlinearity"
    )
  )
}

# ---------------------------------------------------------
# 4. Define input paths
# ---------------------------------------------------------

file_input_rds <- file.path(dir_data_derived, "ms4_analysis_derived.rds")
file_input_csv <- file.path(dir_data_derived, "ms4_analysis_derived.csv")

file_exp1_dataset <- file.path(dir_tables, "06_exp1_light_dataset_summary.csv")
file_exp2_dataset <- file.path(dir_tables, "07_exp2_defined_particles_dataset_summary.csv")
file_exp3_dataset <- file.path(dir_tables, "08_exp3_brake_size_dataset_summary.csv")
file_exp4_dataset <- file.path(dir_tables, "09_exp4_field_spm_dataset_summary.csv")

file_exp1_summary <- file.path(dir_tables, "06_exp1_light_model_summary_table.csv")
file_exp2_summary <- file.path(dir_tables, "07_exp2_defined_particles_model_summary_table.csv")
file_exp3_summary <- file.path(dir_tables, "08_exp3_brake_size_model_summary_table.csv")
file_exp4_summary <- file.path(dir_tables, "09_exp4_field_spm_model_summary_table.csv")

file_exp1_compare <- file.path(dir_tables, "06_exp1_light_model_comparison.csv")
file_exp2_compare <- file.path(dir_tables, "07_exp2_defined_particles_model_comparison.csv")
file_exp3_compare <- file.path(dir_tables, "08_exp3_brake_size_model_comparison.csv")
file_exp4_compare <- file.path(dir_tables, "09_exp4_field_spm_model_comparison.csv")

file_exp4_day4_compare <- file.path(
  project_root, "outputs", "models", "exp4_field_spm",
  "exp4_field_spm_day4_model_comparison.csv"
)

timestamp_now <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_log <- file.path(dir_logs, paste0("11_tables_main_log_", timestamp_now, ".txt"))
file_manifest <- file.path(dir_tables, "11_tables_main_manifest.csv")

# ---------------------------------------------------------
# 5. Load derived dataset
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
      file_input_csv
    ),
    call. = FALSE
  )
}

cat("Derived dataset loaded from:\n")
cat(input_source_used, "\n\n")

# ---------------------------------------------------------
# 6. Load saved outputs
# ---------------------------------------------------------

exp1_dataset <- read_if_exists(file_exp1_dataset)
exp2_dataset <- read_if_exists(file_exp2_dataset)
exp3_dataset <- read_if_exists(file_exp3_dataset)
exp4_dataset <- read_if_exists(file_exp4_dataset)

exp1_summary <- read_if_exists(file_exp1_summary)
exp2_summary <- read_if_exists(file_exp2_summary)
exp3_summary <- read_if_exists(file_exp3_summary)
exp4_summary <- read_if_exists(file_exp4_summary)

exp1_compare <- read_if_exists(file_exp1_compare)
exp2_compare <- read_if_exists(file_exp2_compare)
exp3_compare <- read_if_exists(file_exp3_compare)
exp4_compare <- read_if_exists(file_exp4_compare)
exp4_day4_compare <- read_if_exists(file_exp4_day4_compare)

manifest_tables <- list()

# ---------------------------------------------------------
# 7. Table 1: experimental design summary
# ---------------------------------------------------------

table1 <- tibble::tribble(
  ~experiment_num, ~experiment_name, ~hypothesis, ~primary_exposure, ~days_sampled, ~design_summary,
  9.2, "Experiment 1: Light-only irradiance gradient", "H1", "lux_exposure", "4, 8, 12, 16", "Particle-free optical gradient to test whether reduced irradiance alone alters motility.",
  10.2, "Experiment 2: Defined particle concentration series", "H2", "log10(NTU + 1)", "4, 8, 12, 16", "Defined mineral and organic particles across controlled NTU treatments to isolate compositional effects.",
  11.2, "Experiment 3: Brake-wear size comparison", "H3", "log10(brake-wear concentration + 1)", "4 only", "Brake-wear coarse versus fine fractions compared under matched mass loading.",
  8.2, "Experiment 4: Field-derived SPM gradient", "H4", "log10(NTU + 1)", "4, 8, 12, 16", "Environmentally realistic whole-mixture field SPM gradient with nonlinear NTU-response explicitly evaluated."
)

manifest_tables[[length(manifest_tables) + 1]] <- save_table_multi(
  data = table1,
  file_stem = "Table1_experimental_design_summary",
  caption_text = paste(
    "Table 1. Experimental design summary for Manuscript 4.",
    "The four experiments were designed to isolate optical, compositional, size-dependent, and environmentally realistic whole-mixture pathways affecting zoospore motility."
  ),
  width = 12
)

cat("Table 1 saved.\n")

# ---------------------------------------------------------
# 8. Table 2: dataset summary
# ---------------------------------------------------------

table2 <- dplyr::bind_rows(
  if (!is.null(exp1_dataset)) dplyr::mutate(exp1_dataset, experiment = "Experiment 1"),
  if (!is.null(exp2_dataset)) dplyr::mutate(exp2_dataset, experiment = "Experiment 2"),
  if (!is.null(exp3_dataset)) dplyr::mutate(exp3_dataset, experiment = "Experiment 3"),
  if (!is.null(exp4_dataset)) dplyr::mutate(exp4_dataset, experiment = "Experiment 4")
) |>
  dplyr::relocate(experiment)

manifest_tables[[length(manifest_tables) + 1]] <- save_table_multi(
  data = table2,
  file_stem = "Table2_dataset_summary",
  caption_text = paste(
    "Table 2. Dataset summary for Manuscript 4 experiments.",
    "Values are drawn from experiment-specific analysis-ready subsets and summarise exposure ranges and response counts used in modelling."
  ),
  width = 13
)

cat("Table 2 saved.\n")

# ---------------------------------------------------------
# 9. Table 3: preferred model summary
# ---------------------------------------------------------

table3 <- dplyr::bind_rows(
  if (!is.null(exp1_summary)) dplyr::mutate(exp1_summary, experiment = "Experiment 1"),
  if (!is.null(exp2_summary)) dplyr::mutate(exp2_summary, experiment = "Experiment 2"),
  if (!is.null(exp3_summary)) dplyr::mutate(exp3_summary, experiment = "Experiment 3"),
  if (!is.null(exp4_summary)) dplyr::mutate(exp4_summary, experiment = "Experiment 4")
) |>
  dplyr::relocate(experiment)

manifest_tables[[length(manifest_tables) + 1]] <- save_table_multi(
  data = table3,
  file_stem = "Table3_model_summary_table",
  caption_text = paste(
    "Table 3. Preferred model summaries for Manuscript 4 experiments.",
    "Fixed-effect estimates are shown for the preferred model selected within each experiment."
  ),
  width = 14
)

cat("Table 3 saved.\n")

# ---------------------------------------------------------
# 10. Table 4: model comparison / sensitivity summary
# ---------------------------------------------------------

table4_general <- dplyr::bind_rows(
  if (!is.null(exp1_compare)) {
    exp1_compare |>
      dplyr::slice(1) |>
      dplyr::transmute(
        experiment = "Experiment 1",
        preferred_model = model,
        preferred_aic = AIC,
        next_best_delta_aic = dplyr::lead(delta_aic, default = NA_real_)[1],
        interpretation = "Evaluate whether light-only irradiance explains motility response."
      )
  },
  if (!is.null(exp2_compare)) {
    exp2_compare |>
      dplyr::slice(1) |>
      dplyr::transmute(
        experiment = "Experiment 2",
        preferred_model = model,
        preferred_aic = AIC,
        next_best_delta_aic = dplyr::lead(delta_aic, default = NA_real_)[1],
        interpretation = "Evaluate whether particle type modifies NTU-response."
      )
  },
  if (!is.null(exp3_compare)) {
    exp3_compare |>
      dplyr::slice(1) |>
      dplyr::transmute(
        experiment = "Experiment 3",
        preferred_model = model,
        preferred_aic = AIC,
        next_best_delta_aic = dplyr::lead(delta_aic, default = NA_real_)[1],
        interpretation = "Evaluate whether brake-wear size modifies concentration-response."
      )
  },
  if (!is.null(exp4_compare)) {
    exp4_compare |>
      dplyr::slice(1) |>
      dplyr::transmute(
        experiment = "Experiment 4",
        preferred_model = model,
        preferred_aic = AIC,
        next_best_delta_aic = dplyr::lead(delta_aic, default = NA_real_)[1],
        interpretation = "Evaluate whether field-derived SPM follows a simple monotonic NTU-response."
      )
  }
)

table4_h4 <- dplyr::bind_rows(
  extract_h4(exp4_compare, "Experiment 4 full dataset"),
  extract_h4(exp4_day4_compare, "Experiment 4 Day 4 sensitivity")
)

table4 <- dplyr::bind_rows(
  table4_general,
  table4_h4
)

manifest_tables[[length(manifest_tables) + 1]] <- save_table_multi(
  data = table4,
  file_stem = "Table4_model_comparison_summary",
  caption_text = paste(
    "Table 4. Model comparison and sensitivity summary for Manuscript 4.",
    "For Experiments 1-4, the preferred model and separation from the next-best candidate are summarised.",
    "For Experiment 4, the nonlinear versus monotonic NTU comparison is shown for both the full dataset and the Day 4 sensitivity analysis."
  ),
  width = 13
)

cat("Table 4 saved.\n")

# ---------------------------------------------------------
# 11. Save manifest and log
# ---------------------------------------------------------

manifest_tables_df <- dplyr::bind_rows(manifest_tables)
readr::write_csv(manifest_tables_df, file_manifest)

sink(file_log)
cat("SPM Analysis - 11_tables_main log\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Project:\n")
cat(project_title, "\n")
cat("Manuscript:\n")
cat(manuscript_short, "\n\n")
cat("Input dataset:\n")
cat(input_source_used, "\n\n")
cat("Tables saved:\n")
print(manifest_tables_df)
cat("\n")
cat("Session information:\n\n")
print(utils::sessionInfo())
sink()

cat("Manifest written to:\n")
cat(file_manifest, "\n")
cat("Log written to:\n")
cat(file_log, "\n\n")

cat("========================================================\n")
cat("SCRIPT 11 COMPLETE: TABLES MAIN\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")