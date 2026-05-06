# =========================================================
# Script title: 11_tables_main-v2.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Generate the main publication tables for
#          Manuscript 4 from the derived dataset and saved
#          experiment-specific model outputs.
# Focus: hypothesis-driven outputs (H1–H4)
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
# Outputs:
# - outputs/tables/Table1_experimental_design_summary.csv
# - outputs/tables/Table2_dataset_summary.csv
# - outputs/tables/Table3_model_summary_table.csv
# - outputs/tables/Table4_model_comparison_summary.csv
# - outputs/tables/11_tables_main_manifest.csv
# - outputs/logs/11_tables_main_log_*.txt
# Date created: 26 March 2026
# Last updated: 26 March 2026
# Notes/dependencies:
# - Run 01_setup_packages_and_paths.R first.
# - Run 04_derive_variables.R before this script.
# - Run 06-09 model scripts before this script.
# - Tables follow the publication plan:
#   Table 1 Experimental design summary
#   Table 2 Dataset summary
#   Table 3 Model summary table
#   Table 4 Model comparison / sensitivity summary
# - The pipeline is manuscript-driven and restricted to the
#   Manuscript 4 experimental block:
#   experiment_num %in% c(8.2, 9.2, 10.2, 11.2)
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 11: TABLES MAIN\n")
cat("========================================================\n\n")

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# --- Load model comparisons ---
mod_exp4 <- read_if_exists("outputs/tables/09_exp4_field_spm_model_comparison.csv")
mod_exp4_day4 <- read_if_exists("outputs/models/exp4_field_spm/exp4_field_spm_day4_model_comparison.csv")

# --- Extract H4 result ---
extract_h4 <- function(df, label) {
  if (is.null(df)) return(NULL)
  
  mono <- df |> dplyr::filter(model == "monotonic_primary") |> dplyr::pull(AIC)
  nonlin <- df |> dplyr::filter(model == "nonlinear_primary") |> dplyr::pull(AIC)
  
  tibble::tibble(
    dataset = label,
    delta_aic = mono - nonlin,
    interpretation = dplyr::case_when(
      (mono - nonlin) > 10 ~ "Strong support for nonlinear NTU-response",
      (mono - nonlin) > 2 ~ "Moderate support",
      TRUE ~ "No strong support"
    )
  )
}

table4 <- dplyr::bind_rows(
  extract_h4(mod_exp4, "Full dataset"),
  extract_h4(mod_exp4_day4, "Day 4 sensitivity")
)

readr::write_csv(table4, "outputs/tables/Table4_H4_summary.csv")

cat("Table 4 (H4 summary) saved.\n")