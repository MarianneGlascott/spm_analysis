# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      02_qc_report.R
# Author:      Marianne Glascott
# Date:        2026-02-10
# Purpose:
#   QC reporting only: generate QC flags + summaries.
#   NO rows are dropped and NO exclusions are applied in this script.
#
# Input:
#   data_processed/ms4_clean.parquet
#
# Outputs (written to disk):
#   outputs/qc_flags_table.csv
#   outputs/qc_summary.txt
#   outputs/qc_missingness_by_variable.csv
#   outputs/qc_missingness_by_group.csv
#   outputs/qc_total_vs_components_mismatches.csv   (should be empty)
#   outputs/qc_proposed_exclusions.csv              (template only; not applied)
#   outputs/02_qc_report_log.txt

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(arrow)
  library(tidyr)
  library(tibble)
})

# ---- Logging ----
log_file <- here::here("outputs", "02_qc_report_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 02_qc_report.R")

# ---- Paths ----
in_path <- here::here("data_processed", "ms4_clean.parquet")

out_flags <- here::here("outputs", "qc_flags_table.csv")
out_sum   <- here::here("outputs", "qc_summary.txt")
out_missv <- here::here("outputs", "qc_missingness_by_variable.csv")
out_missg <- here::here("outputs", "qc_missingness_by_group.csv")
out_misct <- here::here("outputs", "qc_total_vs_components_mismatches.csv")
out_excl  <- here::here("outputs", "qc_proposed_exclusions.csv")

if (!file.exists(in_path)) {
  stop("Input parquet not found at: ", in_path, "\nRun scripts/01_import_clean.R first.")
}

df <- arrow::read_parquet(in_path)
log_msg(paste0("Loaded parquet: n=", nrow(df), " rows; p=", ncol(df), " columns."))

# ---- Confirmed QC decisions (from Marianne) ----
# QC-1: total_cells should equal mobile_cell_count + stationary_cell_count
#       We treat mismatch as a QC flag for inspection (not automatic exclusion).
# QC-2: missing experiment is allowed (NOT a QC issue).
# QC-3: missing tile_count is allowed (NOT a QC issue).

log_msg("QC decisions:")
log_msg(" - QC-1: check total == mobile+stationary (flag + save mismatches; no exclusions).")
log_msg(" - QC-2: missing experiment allowed (no flag).")
log_msg(" - QC-3: missing tile_count allowed (no flag).")

# ---- QC-1 diagnostic table: total vs components ----
qc_total_mismatch <- df %>%
  filter(
    !is.na(total_cells),
    !is.na(mobile_cell_count),
    !is.na(stationary_cell_count),
    total_cells != mobile_cell_count + stationary_cell_count
  ) %>%
  mutate(diff = total_cells - (mobile_cell_count + stationary_cell_count)) %>%
  select(video_file, species, site, season, culture, experiment, well,
         total_cells, mobile_cell_count, stationary_cell_count, diff) %>%
  arrange(desc(abs(diff)))

# Always write QC-1 mismatch file, even if empty
if (nrow(qc_total_mismatch) == 0) {
  readr::write_csv(
    tibble::tibble(
      video_file = character(),
      species = character(),
      site = character(),
      season = character(),
      culture = character(),
      experiment = character(),
      well = character(),
      total_cells = integer(),
      mobile_cell_count = integer(),
      stationary_cell_count = integer(),
      diff = integer()
    ),
    out_misct
  )
  log_msg("QC-1 mismatch table written (0 rows; no mismatches detected).")
} else {
  readr::write_csv(qc_total_mismatch, out_misct)
  log_msg(paste0("QC-1 mismatch table written (n=", nrow(qc_total_mismatch), ")."))
}


# ---- Row-level QC flags (no exclusion) ----
flags <- df %>%
  mutate(
    qc_flag_total_cells_missing = is.na(total_cells),
    qc_flag_total_cells_nonpositive = !is.na(total_cells) & total_cells <= 0,
    
    qc_flag_mobile_missing = is.na(mobile_cell_count),
    qc_flag_stationary_missing = is.na(stationary_cell_count),
    
    # QC-1: mismatch flag (should be FALSE for all rows currently)
    qc_flag_total_mismatch = ifelse(
      !is.na(total_cells) & !is.na(mobile_cell_count) & !is.na(stationary_cell_count),
      total_cells != mobile_cell_count + stationary_cell_count,
      FALSE
    ),
    
    qc_flag_non_mobile_negative = !is.na(non_mobile_count) & non_mobile_count < 0,
    
    qc_flag_motility_ratio_missing = is.na(motility_ratio),
    qc_flag_motility_ratio_outside_0_1 = !is.na(motility_ratio) & (motility_ratio < 0 | motility_ratio > 1),
    
    qc_flag_days_from_start_missing = is.na(days_from_start),
    qc_flag_ntu_missing = is.na(ntu)
    
    # NOTE: no flags for missing experiment or tile_count (explicitly allowed)
  )

flag_cols <- names(flags)[grepl("^qc_flag_", names(flags))]
flags <- flags %>%
  mutate(qc_any_flag = Reduce(`|`, across(all_of(flag_cols)), init = FALSE))

# ---- Write QC flags table ----
id_cols <- intersect(
  c("video_file","camera_file","species","site","season","culture","experiment","well",
    "count_date_yyyy_mm_dd","start_date_yyyy_mm_dd","days_from_start"),
  names(flags)
)

flags_out <- flags %>%
  select(all_of(id_cols), all_of(flag_cols), qc_any_flag)

readr::write_csv(flags_out, out_flags)
log_msg(paste0("Wrote qc flags table: ", out_flags))

# ---- Missingness summaries ----
missing_by_var <- tibble(
  variable = names(df),
  n_missing = vapply(df, function(x) sum(is.na(x)), integer(1)),
  prop_missing = n_missing / nrow(df)
) %>%
  arrange(desc(n_missing))

readr::write_csv(missing_by_var, out_missv)
log_msg(paste0("Wrote missingness by variable: ", out_missv))

group_vars <- intersect(c("site","species","season"), names(df))
if (length(group_vars) > 0) {
  missing_by_group <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = n(),
      ntu_missing = sum(is.na(ntu)),
      days_from_start_missing = sum(is.na(days_from_start)),
      total_mismatch = sum(qc_flag_total_mismatch),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  readr::write_csv(missing_by_group, out_missg)
  log_msg(paste0("Wrote missingness by group: ", out_missg))
} else {
  log_msg("Skipped missingness-by-group (no grouping vars found).")
}

# ---- QC summary text ----
flag_counts <- tibble(
  flag = flag_cols,
  n_true = vapply(flag_cols, function(fc) sum(flags[[fc]] %in% TRUE, na.rm = TRUE), integer(1))
) %>%
  arrange(desc(n_true))

summary_lines <- c(
  "Manuscript 4 — 02_qc_report.R summary",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Input: ", in_path),
  paste0("Rows: ", nrow(df)),
  paste0("Columns: ", ncol(df)),
  "",
  "QC decisions (explicit):",
  "- QC-1: check total == mobile+stationary (mismatches saved; no exclusions)",
  "- QC-2: missing experiment allowed (no QC flag)",
  "- QC-3: missing tile_count allowed (no QC flag)",
  "",
  paste0("QC-1 mismatches found: ", nrow(qc_total_mismatch)),
  "",
  "QC flag counts (TRUE):"
)

summary_lines <- c(summary_lines, capture.output(print(flag_counts, n = Inf, width = Inf)))
summary_lines <- c(summary_lines, "", paste0("Rows with ANY qc flag TRUE: ", sum(flags$qc_any_flag, na.rm = TRUE)))

writeLines(summary_lines, out_sum)
log_msg(paste0("Wrote QC summary: ", out_sum))

# ---- Proposed exclusions template (NOT applied) ----
excl_template <- tibble(
  video_file = character(),
  reason = character(),
  rule_id = character(),
  proposed_by = character(),
  proposed_on = character()
)
readr::write_csv(excl_template, out_excl)
log_msg(paste0("Wrote exclusions template (empty; not applied): ", out_excl))

log_msg("02_qc_report.R complete (no rows excluded).")


