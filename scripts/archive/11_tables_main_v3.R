# =========================================================
# Script title: 11_tables_main_v3.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Generate the main publication tables for
#          Manuscript 4 from the derived dataset and saved
#          experiment-specific model outputs.
# Focus: hypothesis-driven outputs (H1–H4)
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
# 2. Small local helpers
# ---------------------------------------------------------

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

write_caption_md_local <- function(table_name, caption_text, subdir = NULL) {
  out_dir <- if (is.null(subdir)) dir_tables else file.path(dir_tables, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  caption_file <- file.path(out_dir, paste0(table_name, "_caption.md"))
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
                             width = 10,
                             height = NULL,
                             base_size = 10) {
  if (!requireNamespace("gridExtra", quietly = TRUE) ||
      !requireNamespace("grid", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Packages 'gridExtra', 'grid', and 'ggplot2' are required for table rendering.",
      call. = FALSE
    )
  }
  
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
  
  # PDF
  grDevices::cairo_pdf(file_pdf, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(tg)
  grDevices::dev.off()
  
  # PNG
  ggplot2::ggsave(
    filename = file_png,
    plot = tg,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  # TIFF
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
    caption_file <- write_caption_md_local(
      table_name = file_stem,
      caption_text = caption_text
    )
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
# 3. Define input paths
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
file_log <- file.path(
  dir_logs,
  paste0("11_tables_main_log_", timestamp_now, ".txt")
)
file_manifest <- file.path(dir_tables, "11_tables_main_manifest.csv")

# ---------------------------------------------------------
# 4. Load main derived dataset
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
# 5. Load saved experiment outputs
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

# ---------------------------------------------------------
# 6. Build Table 1: experimental design summary
# ---------------------------------------------------------

table1 <- tibble::tribble(
  ~experiment_num, ~experiment_name, ~hypothesis, ~primary_exposure, ~design_summary,
  9.2, "Experiment 1: Light-only irradiance gradient", "H1", "Lux / irradiance", "Particle-free optical gradient to test whether reduced light alone explains motility response.",
  10.2, "Experiment 2: Defined particle concentration", "H2", "log10(NTU + 1)", "Defined mineral and organic particles across controlled NTU treatments to isolate compositional effects.",
  11.2, "Experiment 3: Brake-wear size under equal mass loading", "H3", "log10(brake-wear concentration + 1)", "Brake-wear coarse versus fine fractions compared under matched mass loading.",
  8.2, "Experiment 4: Field-derived SPM gradient", "H4", "log10(NTU + 1)", "Environmentally realistic whole-mixture field SPM gradient with nonlinear NTU-response explicitly evaluated."
)

manifest_parts <- list()

manifest_parts[[length(manifest_parts) + 1]] <- save_table_multi(
  data = table1,
  file_stem = "Table1_experimental_design_summary",
  caption_text = paste(
    "Table 1. Experimental design summary for Manuscript 4.",
    "The four experiments were designed to isolate optical, compositional, size-dependent, and environmentally realistic whole-mixture pathways affecting zoospore motility."
  ),
  width = 11
)

cat("Table 1 saved.\n")

# ---------------------------------------------------------
# 7. Build Table 2: dataset summary
# ---------------------------------------------------------

table2 <- dplyr::bind_rows(
  if (!is.null(exp1_dataset)) dplyr::mutate(exp1_dataset, experiment = "Experiment 1"),
  if (!is.null(exp2_dataset)) dplyr::mutate(exp2_dataset, experiment = "Experiment 2"),
  if (!is.null(exp3_dataset)) dplyr::mutate(exp3_dataset, experiment = "Experiment 3"),
  if (!is.null(exp4_dataset)) dplyr::mutate(exp4_dataset, experiment = "Experiment 4")
) |>
  dplyr::relocate(experiment)

manifest_parts[[length(manifest_parts) + 1]] <- save_table_multi(
  data = table2,
  file_stem = "Table2_dataset_summary",
  caption_text = paste(
    "Table 2. Dataset summary for Manuscript 4 experiments.",
    "Values are drawn from experiment-specific analysis-ready subsets and summarise the observed exposure ranges and response counts used in modelling."
  ),
  width = 12
)

cat("Table 2 saved.\n")

# ---------------------------------------------------------
# 8. Build Table 3: preferred model summary table
# ---------------------------------------------------------

table3 <- dplyr::bind_rows(
  if (!is.null(exp1_summary)) dplyr::mutate(exp1_summary, experiment = "Experiment 1"),
  if (!is.null(exp2_summary)) dplyr::mutate(exp2_summary, experiment = "Experiment 2"),
  if (!is.null(exp3_summary)) dplyr::mutate(exp3_summary, experiment = "Experiment 3"),
  if (!is.null(exp4_summary)) dplyr::mutate(exp4_summary, experiment = "Experiment 4")
) |>
  dplyr::relocate(experiment)

manifest_parts[[length(manifest_parts) + 1]] <- save_table_multi(
  data = table3,
  file_stem = "Table3_model_summary_table",
  caption_text = paste(
    "Table 3. Preferred model summaries for Manuscript 4 experiments.",
    "Fixed-effect estimates are shown for the preferred model selected within each experiment."
  ),
  width = 13
)

cat("Table 3 saved.\n")

# ---------------------------------------------------------
# 9. Build Table 4: model comparison / sensitivity summary
# ---------------------------------------------------------

table4_core <- dplyr::bind_rows(
  extract_h4(exp4_compare, "Full dataset"),
  extract_h4(exp4_day4_compare, "Day 4 sensitivity")
)

table4_other <- dplyr::bind_rows(
  if (!is.null(exp1_compare)) dplyr::mutate(exp1_compare, experiment = "Experiment 1"),
  if (!is.null(exp2_compare)) dplyr::mutate(exp2_compare, experiment = "Experiment 2"),
  if (!is.null(exp3_compare)) dplyr::mutate(exp3_compare, experiment = "Experiment 3"),
  if (!is.null(exp4_compare)) dplyr::mutate(exp4_compare, experiment = "Experiment 4")
)

manifest_parts[[length(manifest_parts) + 1]] <- save_table_multi(
  data = table4_core,
  file_stem = "Table4_model_comparison_summary",
  caption_text = paste(
    "Table 4. Experiment 4 model comparison and sensitivity summary.",
    "The table compares monotonic and nonlinear NTU-response forms for the full field-derived SPM dataset and the Day 4 sensitivity analysis."
  ),
  width = 10
)

# Optional: also save the broader model comparison table for traceability
manifest_parts[[length(manifest_parts) + 1]] <- save_table_multi(
  data = table4_other,
  file_stem = "Table4a_all_model_comparisons",
  caption_text = paste(
    "Table 4a. Model comparison outputs across experiments.",
    "AIC-based comparisons are shown for the candidate models fitted within each experiment."
  ),
  width = 12
)

cat("Table 4 saved.\n")

# ---------------------------------------------------------
# 10. Write manifest
# ---------------------------------------------------------

manifest <- dplyr::bind_rows(manifest_parts)
readr::write_csv(manifest, file_manifest)

cat("Manifest written to:\n")
cat(file_manifest, "\n\n")

# ---------------------------------------------------------
# 11. Write log
# ---------------------------------------------------------

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
print(manifest)
cat("\n")

cat("Session information:\n\n")
print(utils::sessionInfo())
sink()

cat("Log written to:\n")
cat(file_log, "\n\n")

cat("========================================================\n")
cat("SCRIPT 11 COMPLETE: TABLES MAIN\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")