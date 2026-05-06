# =========================================================
# Script title: 00_run_all.R
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Project: SPM Analysis
# Manuscript: Manuscript 4
# Purpose: Optional master script to execute the full
#          analysis pipeline end-to-end in numbered order.
# Inputs: Project directory structure and numbered scripts
#         in /scripts
# Outputs: All standard derived data, models, figures,
#          tables, logs, and session information written
#          to disk by downstream scripts
# Date created: 24 March 2026
# Last updated: 24 March 2026
# Notes/dependencies:
# - This script is optional.
# - It is intended for full reproducibility, archiving,
#   and end-to-end reruns.
# - During active development, scripts may still be run
#   individually in sequence.
# - This script contains orchestration only and does not
#   contain substantive analysis code.
# =========================================================

rm(list = ls())
graphics.off()

cat("\n========================================================\n")
cat("SPM ANALYSIS: FULL PIPELINE START\n")
cat("Manuscript 4 - Kelp zoospore motility\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check required packages for bootstrapping the project
# ---------------------------------------------------------

required_bootstrap_packages <- c("here", "renv")

missing_bootstrap_packages <- required_bootstrap_packages[
  !vapply(required_bootstrap_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bootstrap_packages) > 0) {
  stop(
    paste0(
      "Missing required bootstrap package(s): ",
      paste(missing_bootstrap_packages, collapse = ", "),
      ". Please install these before running 00_run_all.R."
    ),
    call. = FALSE
  )
}

# Activate project-local package environment where available
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
  cat("renv activated successfully.\n")
} else {
  cat("No renv/activate.R found. Proceeding without explicit renv activation.\n")
}

# ---------------------------------------------------------
# 2. Confirm project root and required directory structure
# ---------------------------------------------------------

project_root <- here::here()

cat("Project root detected as:\n")
cat(project_root, "\n\n")

required_dirs <- c(
  "data_raw",
  "data_clean",
  "data_derived",
  "scripts",
  "R",
  "outputs",
  "outputs/figures",
  "outputs/tables",
  "outputs/models",
  "outputs/logs",
  "reports"
)

missing_dirs <- required_dirs[!dir.exists(file.path(project_root, required_dirs))]

if (length(missing_dirs) > 0) {
  stop(
    paste0(
      "The following required directory/directories are missing:\n- ",
      paste(missing_dirs, collapse = "\n- "),
      "\nPlease create them before running the full pipeline."
    ),
    call. = FALSE
  )
}

cat("Required project directories verified.\n\n")

# ---------------------------------------------------------
# 3. Define ordered script pipeline
# ---------------------------------------------------------

script_paths <- c(
  "scripts/01_setup_packages_and_paths.R",
  "scripts/02_import_and_clean_data.R",
  "scripts/03_qc_and_exclusions.R",
  "scripts/04_derive_variables.R",
  "scripts/05_eda.R",
  "scripts/06_models_exp1_light.R",
  "scripts/07_models_exp2_defined_particles.R",
  "scripts/08_models_exp3_brake_size.R",
  "scripts/09_models_exp4_field_spm.R",
  "scripts/10_figures_main.R",
  "scripts/11_tables_main.R",
  "scripts/12_supplementary_outputs.R"
)

missing_scripts <- script_paths[!file.exists(file.path(project_root, script_paths))]

if (length(missing_scripts) > 0) {
  stop(
    paste0(
      "The following pipeline script(s) are missing:\n- ",
      paste(missing_scripts, collapse = "\n- "),
      "\nPlease check the scripts/ folder before running the pipeline."
    ),
    call. = FALSE
  )
}

cat("All pipeline scripts found.\n\n")

# ---------------------------------------------------------
# 4. Helper object for runtime logging
# ---------------------------------------------------------

pipeline_log <- data.frame(
  step_number = integer(),
  script = character(),
  start_time = character(),
  end_time = character(),
  elapsed_seconds = numeric(),
  status = character(),
  stringsAsFactors = FALSE
)

log_file_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
pipeline_log_path <- file.path(
  project_root,
  "outputs",
  "logs",
  paste0("pipeline_run_log_", log_file_timestamp, ".csv")
)

session_info_path <- file.path(
  project_root,
  "outputs",
  "logs",
  paste0("session_info_", log_file_timestamp, ".txt")
)

# ---------------------------------------------------------
# 5. Run scripts in sequence
# ---------------------------------------------------------

for (i in seq_along(script_paths)) {
  current_script <- script_paths[i]
  script_start <- Sys.time()

  cat("--------------------------------------------------------\n")
  cat("Running step", i, "of", length(script_paths), "\n")
  cat("Script:", current_script, "\n")
  cat("Started:", format(script_start, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("--------------------------------------------------------\n")

  result <- tryCatch(
    {
      source(file.path(project_root, current_script), local = FALSE, echo = FALSE)
      list(status = "success", message = NA_character_)
    },
    error = function(e) {
      list(status = "error", message = conditionMessage(e))
    }
  )

  script_end <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(script_end, script_start, units = "secs"))

  pipeline_log <- rbind(
    pipeline_log,
    data.frame(
      step_number = i,
      script = current_script,
      start_time = format(script_start, "%Y-%m-%d %H:%M:%S"),
      end_time = format(script_end, "%Y-%m-%d %H:%M:%S"),
      elapsed_seconds = round(elapsed_seconds, 2),
      status = result$status,
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(pipeline_log, pipeline_log_path, row.names = FALSE)

  if (identical(result$status, "success")) {
    cat("Completed successfully in", round(elapsed_seconds, 2), "seconds.\n\n")
  } else {
    cat("ERROR in:", current_script, "\n")
    cat("Message:", result$message, "\n\n")

    cat("Pipeline stopped.\n")
    cat("A partial run log has been written to:\n")
    cat(pipeline_log_path, "\n\n")

    stop(
      paste0(
        "Pipeline terminated because ",
        current_script,
        " failed with error: ",
        result$message
      ),
      call. = FALSE
    )
  }
}

# ---------------------------------------------------------
# 6. Save session information for reproducibility
# ---------------------------------------------------------

cat("Writing session information...\n")

sink(session_info_path)
cat("SPM Analysis - session information\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(utils::sessionInfo())
sink()

# ---------------------------------------------------------
# 7. Final completion message
# ---------------------------------------------------------

pipeline_end_time <- Sys.time()
total_elapsed <- sum(pipeline_log$elapsed_seconds, na.rm = TRUE)

cat("\n========================================================\n")
cat("SPM ANALYSIS: FULL PIPELINE COMPLETE\n")
cat("End time:", format(pipeline_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total elapsed time (seconds):", round(total_elapsed, 2), "\n")
cat("Run log written to:\n")
cat(pipeline_log_path, "\n")
cat("Session info written to:\n")
cat(session_info_path, "\n")
cat("========================================================\n\n")