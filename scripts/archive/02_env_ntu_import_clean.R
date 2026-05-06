# ==============================================================================
# Project:     SPM Analysis — Environmental Turbidity Context (Paper 4/5)
# Script:      02_env_ntu_import_clean.R
# Author:      Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2026-02-16
# Version:     0.1
#
# Description:
#   Imports and consolidates environmental turbidity logger data sets from multiple
#   CSV files into a single, analysis-ready table with standardized column names,
#   time_stamps, and core metadata.
#
# Instrument families supported:
#   1) Aquatec turbidity logger (NTU measured directly)
#   2) PEBL-CIC Grow Probes (IR proxy; NTU conversion done later)
#
# Rules:
#   - No silent filtering: any exclusions must be logged; prefer flagging over dropping
#   - Preserve raw values + provenance: keep source_file, instrument_type, instrument_id
#   - Standardize time to POSIXct; keep timezone used
#
# Inputs:
#   data_raw/turbidity/aquatec/        (0+ CSV files)
#   data_raw/turbidity/grow_probe/     (0+ CSV files)
#   data_raw/turbidity/metadata/site_registry.csv   
#   data_raw/turbidity/metadata/instrument_registry.csv
#
# Outputs:
#   data_processed/env_ntu_merged.parquet
#   data_processed/env_ntu_merged.csv (optional)
#   outputs/env_ntu_import_log.txt
#   outputs/env_ntu_file_manifest.csv
#   outputs/env_ntu_column_map.csv
#   outputs/env_ntu_quick_summary_by_file.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(janitor)
  library(purrr)
  library(tidyr)
  library(arrow)
})

# ---- 1) Logging ----
log_file <- here::here("outputs", "env_ntu_import_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 02_env_ntu_import_clean.R")

# ---- 2) Paths ----
dir.create(here::here("data_processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs"), recursive = TRUE, showWarnings = FALSE)

aquatec_dir <- here::here("data_raw", "turbidity", "aquatec")
grow_dir    <- here::here("data_raw", "turbidity", "grow_probe")
meta_dir    <- here::here("data_raw", "turbidity", "metadata")

out_parquet <- here::here("data_processed", "env_ntu_merged.parquet")
out_csv     <- here::here("data_processed", "env_ntu_merged.csv")
out_manifest <- here::here("outputs", "env_ntu_file_manifest.csv")
out_colmap   <- here::here("outputs", "env_ntu_column_map.csv")
out_quick    <- here::here("outputs", "env_ntu_quick_summary_by_file.csv")

# ---- 3) Optional metadata (defined outside R) ----
site_registry_path <- here::here(meta_dir, "site_registry.csv")
instr_registry_path <- here::here(meta_dir, "instrument_registry.csv")

site_registry <- NULL
if (file.exists(site_registry_path)) {
  site_registry <- readr::read_csv(site_registry_path, show_col_types = FALSE) %>% clean_names()
  log_msg(paste0("Loaded site registry: ", site_registry_path, " (n=", nrow(site_registry), ")"))
} else {
  log_msg("NOTE: site_registry.csv not found (site fields may remain NA).")
}

instr_registry <- NULL
if (file.exists(instr_registry_path)) {
  instr_registry <- readr::read_csv(instr_registry_path, show_col_types = FALSE) %>% clean_names()
  log_msg(paste0("Loaded instrument registry: ", instr_registry_path, " (n=", nrow(instr_registry), ")"))
} else {
  log_msg("NOTE: instrument_registry.csv not found (instrument_id inference may be weaker).")
}

# ---- 4) Discover files + build manifest ----
list_csv <- function(path) {
  if (!dir.exists(path)) return(character(0))
  list.files(path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
}

aquatec_files <- list_csv(aquatec_dir)
grow_files    <- list_csv(grow_dir)

log_msg(paste0("Discovered Aquatec CSVs: ", length(aquatec_files)))
log_msg(paste0("Discovered Grow Probe CSVs: ", length(grow_files)))

if (length(aquatec_files) + length(grow_files) == 0) {
  stop("No turbidity CSV files found. Check folder paths:\n",
       " - ", aquatec_dir, "\n",
       " - ", grow_dir)
}

manifest <- tibble(
  source_file = c(aquatec_files, grow_files),
  file_name   = basename(source_file),
  instrument_type = c(rep("aquatec", length(aquatec_files)),
                      rep("grow_probe", length(grow_files)))
) %>%
  mutate(
    # basic guesses from filename (improve via registry below)
    instrument_guess = str_extract(file_name, "(GP\\d+|\\d{3,4}[-_]\\d{2,4}|\\d{3,6})"),
    site_guess = str_extract(file_name, "([A-Z]{2,4})") # crude; replace with your own rules if needed
  )

# If instrument registry exists, use patterns to assign instrument_id and maybe site_code
if (!is.null(instr_registry) && all(c("pattern","instrument_id","instrument_type") %in% names(instr_registry))) {
  assign_from_registry <- function(fname, itype) {
    rr <- instr_registry %>%
      filter(instrument_type == itype) %>%
      mutate(hit = str_detect(fname, pattern)) %>%
      filter(hit)
    if (nrow(rr) == 0) return(list(instrument_id = NA_character_, site_code = NA_character_, timezone = NA_character_))
    # If multiple match, take first (we log this situation later)
    list(
      instrument_id = rr$instrument_id[[1]],
      site_code     = if ("site_code" %in% names(rr)) rr$site_code[[1]] else NA_character_,
      timezone      = if ("timezone" %in% names(rr)) rr$timezone[[1]] else NA_character_
    )
  }
  
  reg_assigned <- pmap(manifest[,c("file_name","instrument_type")], assign_from_registry)
  reg_df <- tibble(
    instrument_id = map_chr(reg_assigned, "instrument_id"),
    site_code_reg = map_chr(reg_assigned, "site_code"),
    timezone_reg  = map_chr(reg_assigned, "timezone")
  )
  manifest <- bind_cols(manifest, reg_df)
} else {
  manifest <- manifest %>%
    mutate(instrument_id = NA_character_, site_code_reg = NA_character_, timezone_reg = NA_character_)
}

# Final manifest fields (prefer registry, else guesses)
manifest <- manifest %>%
  mutate(
    instrument_id = coalesce(instrument_id, instrument_guess),
    site_code     = coalesce(site_code_reg, site_guess),
    timezone_used = coalesce(timezone_reg, "Europe/London") # default; change if needed
  ) %>%
  select(-instrument_guess, -site_guess)

readr::write_csv(manifest, out_manifest)
log_msg(paste0("Wrote file manifest: ", out_manifest))

# ---- 5) Time parsing helper (flag, don’t drop) ----
parse_time_multi <- function(x, tz = "Europe/London") {
  if (inherits(x, "POSIXct")) return(x)
  if (is.factor(x)) x <- as.character(x)
  x <- str_trim(x)
  
  # attempt several common formats; keep NA if cannot parse
  parsed <- suppressWarnings(parse_date_time(
    x,
    orders = c(
      "dmy HMS", "dmy HM", "dmy",
      "ymd HMS", "ymd HM", "ymd",
      "mdy HMS", "mdy HM", "mdy",
      "d/m/Y HMS", "d/m/Y HM",
      "Y-m-d H:M:S", "Y-m-d H:M"
    ),
    tz = tz
  ))
  
  parsed
}

# ---- 6) Import functions per instrument type ----
# IMPORTANT: these are “tolerant” parsers: we standardise what we can and
# flag parse failures. We don’t throw away rows silently.

import_aquatec <- function(path, instrument_id, site_code, tz) {
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  raw_names <- names(raw)
  df <- raw %>% clean_names()
  
  # Heuristics: find a likely datetime column
  time_col <- names(df)[str_detect(names(df), "time|date")]
  time_col <- time_col[1] %||% NA_character_
  
  # Heuristics: find NTU column
  ntu_col <- names(df)[str_detect(names(df), "ntu|turbid")]
  ntu_col <- ntu_col[1] %||% NA_character_
  
  # Standardise
  out <- df %>%
    mutate(
      instrument_type = "aquatec",
      instrument_id   = instrument_id,
      site_code       = site_code,
      timezone_used   = tz,
      source_file     = path,
      file_name       = basename(path),
      datetime_raw    = if (!is.na(time_col)) .data[[time_col]] else NA,
      datetime        = parse_time_multi(datetime_raw, tz = tz),
      ntu_raw         = if (!is.na(ntu_col)) suppressWarnings(as.numeric(.data[[ntu_col]])) else NA_real_
    ) %>%
    mutate(
      parse_ok_time = !is.na(datetime),
      parse_ok_ntu  = !is.na(ntu_raw)
    )
  
  # Column map for audit (raw -> used)
  colmap <- tibble(
    file_name = basename(path),
    instrument_type = "aquatec",
    datetime_col_used = time_col,
    ntu_col_used = ntu_col,
    n_raw_cols = length(raw_names)
  )
  
  list(data = out, colmap = colmap)
}

import_grow_probe <- function(path, instrument_id, site_code, tz) {
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  raw_names <- names(raw)
  df <- raw %>% clean_names()
  
  time_col <- names(df)[str_detect(names(df), "^time$|time|date")]
  time_col <- time_col[1] %||% NA_character_
  
  # Grow probes typically: IR / irradiance proxy, lux, temp
  ir_col <- names(df)[str_detect(names(df), "^ir$|irradiance|ir1|ir_av|ir_avg")]
  ir_col <- ir_col[1] %||% NA_character_
  
  lux_col <- names(df)[str_detect(names(df), "lux")]
  lux_col <- lux_col[1] %||% NA_character_
  
  temp_col <- names(df)[str_detect(names(df), "temp")]
  temp_col <- temp_col[1] %||% NA_character_
  
  out <- df %>%
    mutate(
      instrument_type = "grow_probe",
      instrument_id   = instrument_id,
      site_code       = site_code,
      timezone_used   = tz,
      source_file     = path,
      file_name       = basename(path),
      datetime_raw    = if (!is.na(time_col)) .data[[time_col]] else NA,
      datetime        = parse_time_multi(datetime_raw, tz = tz),
      ir_raw          = if (!is.na(ir_col)) suppressWarnings(as.numeric(.data[[ir_col]])) else NA_real_,
      lux_raw         = if (!is.na(lux_col)) suppressWarnings(as.numeric(.data[[lux_col]])) else NA_real_,
      temp_c          = if (!is.na(temp_col)) suppressWarnings(as.numeric(.data[[temp_col]])) else NA_real_
    ) %>%
    mutate(
      parse_ok_time = !is.na(datetime),
      parse_ok_ir   = !is.na(ir_raw)
    )
  
  colmap <- tibble(
    file_name = basename(path),
    instrument_type = "grow_probe",
    datetime_col_used = time_col,
    ir_col_used = ir_col,
    lux_col_used = lux_col,
    temp_col_used = temp_col,
    n_raw_cols = length(raw_names)
  )
  
  list(data = out, colmap = colmap)
}

`%||%` <- function(a, b) if (length(a) == 0 || all(is.na(a))) b else a

# ---- 7) Import all files (with logging) ----
log_msg("Importing Aquatec files...")
aqu_list <- pmap(
  list(manifest$source_file[manifest$instrument_type=="aquatec"],
       manifest$instrument_id[manifest$instrument_type=="aquatec"],
       manifest$site_code[manifest$instrument_type=="aquatec"],
       manifest$timezone_used[manifest$instrument_type=="aquatec"]),
  ~ import_aquatec(..1, ..2, ..3, ..4)
)

log_msg("Importing Grow Probe files...")
gp_list <- pmap(
  list(manifest$source_file[manifest$instrument_type=="grow_probe"],
       manifest$instrument_id[manifest$instrument_type=="grow_probe"],
       manifest$site_code[manifest$instrument_type=="grow_probe"],
       manifest$timezone_used[manifest$instrument_type=="grow_probe"]),
  ~ import_grow_probe(..1, ..2, ..3, ..4)
)

aqu_df <- bind_rows(map(aqu_list, "data"))
gp_df  <- bind_rows(map(gp_list,  "data"))

colmap_df <- bind_rows(
  bind_rows(map(aqu_list, "colmap")),
  bind_rows(map(gp_list,  "colmap"))
)

readr::write_csv(colmap_df, out_colmap)
log_msg(paste0("Wrote column map: ", out_colmap))

# ---- 8) Harmonise schemas (keep superset; don’t drop columns) ----
# Standard set: datetime, instrument_id, site_code, ntu_raw (Aquatec), ir_raw (Grow probe), etc.
merged <- bind_rows(aqu_df, gp_df) %>%
  mutate(
    # optional: join site details if provided
    site_code = na_if(site_code, ""),
    instrument_id = na_if(instrument_id, "")
  )

if (!is.null(site_registry) && "site_code" %in% names(site_registry)) {
  merged <- merged %>% left_join(site_registry, by = "site_code")
  log_msg("Joined site registry onto merged data.")
}

# ---- 9) Minimal integrity checks (flags + summaries) ----
n_total <- nrow(merged)
n_time_bad <- sum(is.na(merged$datetime))
log_msg(paste0("Merged rows: ", n_total))
log_msg(paste0("Rows with unparsed datetime (flagged, not dropped): ", n_time_bad))

# Quick summary by file (for debugging)
quick <- merged %>%
  group_by(instrument_type, instrument_id, file_name) %>%
  summarise(
    n_rows = n(),
    dt_min = suppressWarnings(min(datetime, na.rm = TRUE)),
    dt_max = suppressWarnings(max(datetime, na.rm = TRUE)),
    n_time_na = sum(is.na(datetime)),
    n_ntu_nonmissing = sum(!is.na(ntu_raw)),
    n_ir_nonmissing  = sum(!is.na(ir_raw)),
    .groups = "drop"
  )

readr::write_csv(quick, out_quick)
log_msg(paste0("Wrote quick summary by file: ", out_quick))

# ---- 10) Write outputs ----
arrow::write_parquet(merged, out_parquet)
log_msg(paste0("Saved merged parquet: ", out_parquet))

# Optional CSV convenience (can be large)
readr::write_csv(merged, out_csv)
log_msg(paste0("Saved merged CSV: ", out_csv))

log_msg("02_env_ntu_import_clean.R complete.")
