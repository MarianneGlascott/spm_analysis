# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      01_import_clean.R
# Author:      Marianne Glascott
# Date:        2026-02-10
# Purpose:
#   Import raw data, enforce canonical names + types, create derived variables.
# Input:
#   data_raw/All_data_cleaned.csv
# Outputs (written to disk; required):
#   data_processed/ms4_clean.parquet
#   outputs/data_dictionary.csv
#   outputs/data_summary.txt
#   outputs/01_import_clean_log.txt
# Rules:
#   - No dropped rows
#   - No silent filtering
#   - Log all coercions and warnings
#   - Proportions are derived; counts retained
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(janitor)
  library(lubridate)
  library(arrow)
  library(tibble)
})

# ---- Logging ----
log_file <- here::here("outputs", "01_import_clean_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 01_import_clean.R")

# ---- Paths ----
in_path  <- here::here("data_raw", "All_data_cleaned.csv")
out_parq <- here::here("data_processed", "ms4_clean.parquet")
out_dict <- here::here("outputs", "data_dictionary.csv")
out_sum  <- here::here("outputs", "data_summary.txt")
out_map  <- here::here("outputs", "name_map_01_import_clean.csv")

dir.create(dirname(out_parq), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_dict), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_sum),  recursive = TRUE, showWarnings = FALSE)

# ---- Guardrail: enforce expected input location ----
if (!file.exists(in_path)) {
  stop(
    "Input file not found at: ", in_path, "\n\n",
    "Please place your raw CSV here:\n",
    "  data_raw/All_data_cleaned.csv\n\n",
    "You mentioned a file at C:/SPM_Analysis/data/All_data_cleaned.csv.\n",
    "Move/copy it into data_raw/ (we do not read from arbitrary paths to keep the pipeline reproducible)."
  )
}

log_msg(paste0("Reading raw CSV: ", in_path))

# ---- Import ----
raw_df <- readr::read_csv(in_path, show_col_types = FALSE, progress = FALSE)
log_msg(paste0("Raw data loaded: n=", nrow(raw_df), " rows; p=", ncol(raw_df), " columns."))

orig_names <- names(raw_df)

# ---- Standardise column names to snake_case ----
df <- raw_df %>% janitor::clean_names()
clean_names <- names(df)

# Write name map (audit trail)
name_map <- tibble(original_name = orig_names, clean_name = clean_names)
readr::write_csv(name_map, out_map)
log_msg(paste0("Wrote name map: ", out_map))

# ---- Canonicalise variable names (explicit rename map) ----
# We only rename known variants/synonyms; we do not fabricate missing columns.
rename_map <- c(
  # identifiers / provenance
  "cama_file"        = "camera_file",
  "cama_file_name"   = "camera_file",
  "cama_file_"       = "camera_file",
  
  # counts (observed some truncated names)
  "mobile_ce"        = "mobile_cell_count",
  "mobile_cell"      = "mobile_cell_count",
  "mobile_cells"     = "mobile_cell_count",
  
  "stationary"       = "stationary_cell_count",
  "stationary_"      = "stationary_cell_count",
  "stationary_cells" = "stationary_cell_count",
  
  "total_cells_"     = "total_cells",
  "total_cell"       = "total_cells",
  
  # dates
  "count_dat"        = "count_date_yyyy_mm_dd",
  "count_date"       = "count_date_yyyy_mm_dd",
  "start_date"       = "start_date_yyyy_mm_dd",
  
  # exposures / predictors
  "cu_ug_l_"         = "cu_ug_l",
  "lux_expos"        = "lux_exposure",
  "lux_exposure_"    = "lux_exposure",
  
  # metadata variants
  "collection_year"  = "year",
  "preparatio"       = "preparation",
  "preparation_"     = "preparation"
)

# Apply rename_map only where the source column exists AND target does not already exist
for (src in names(rename_map)) {
  tgt <- rename_map[[src]]
  if (src %in% names(df) && !(tgt %in% names(df))) {
    df <- df %>% dplyr::rename(!!tgt := !!src)
    log_msg(paste0("Renamed column: '", src, "' -> '", tgt, "'"))
  }
}

# ---- Canonical ID: video_file ----
# Prefer video_file if already present; otherwise rename video_file_name -> video_file
if ("video_file_name" %in% names(df)) {
  if (!"video_file" %in% names(df)) {
    df <- dplyr::rename(df, video_file = video_file_name)
    log_msg("Renamed column: 'video_file_name' -> 'video_file' (canonical observation ID).")
  } else {
    log_msg("NOTE: Both 'video_file' and 'video_file_name' exist; keeping 'video_file' as canonical ID.")
  }
}

if ("video_file" %in% names(df)) {
  df$video_file <- factor(df$video_file)
  log_msg("Canonicalised 'video_file' as factor (unique observation ID).")
}

# ---- Expected canonical columns (for reporting only) ----
expected_cols <- c(
  "video_file","camera_file","frame_count",
  "mobile_cell_count","stationary_cell_count","total_cells",
  "motility_ratio","count_date_yyyy_mm_dd","start_date_yyyy_mm_dd",
  "species","site","culture","experiment","well","days_from_start",
  "toxin_exposure","cu_ug_l","ntu","tile","tile_count","lux_exposure",
  "comment","collection_date","year","season","collection_site",
  "preparation","release"
)

missing_expected <- setdiff(expected_cols, names(df))
if (length(missing_expected) > 0) {
  log_msg(paste0(
    "WARNING: Missing expected column(s): ",
    paste(missing_expected, collapse = ", "),
    ". No imputation performed."
  ))
}

# ---- Coercion helpers (log NA introduction) ----
coerce_integer <- function(x, nm) {
  before_na <- sum(is.na(x))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- str_trim(x)
  out <- suppressWarnings(as.integer(x))
  after_na <- sum(is.na(out))
  introduced <- after_na - before_na
  if (introduced > 0) log_msg(paste0("WARNING: integer coercion introduced ", introduced, " new NA(s) in '", nm, "'."))
  out
}

coerce_numeric <- function(x, nm) {
  before_na <- sum(is.na(x))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- str_trim(x)
  out <- suppressWarnings(as.numeric(x))
  after_na <- sum(is.na(out))
  introduced <- after_na - before_na
  if (introduced > 0) log_msg(paste0("WARNING: numeric coercion introduced ", introduced, " new NA(s) in '", nm, "'."))
  out
}

coerce_factor <- function(x, nm) {
  if (is.factor(x)) return(x)
  if (is.numeric(x) || is.integer(x)) return(factor(as.character(x)))
  factor(x)
}

coerce_date_dmy <- function(x, nm) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.factor(x)) x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "none", "None", "NA", "N/A")] <- NA_character_
  
  out <- suppressWarnings(lubridate::dmy(x, quiet = TRUE))
  fail <- sum(!is.na(x) & is.na(out))
  if (fail > 0) log_msg(paste0("WARNING: date parsing failed for ", fail, " value(s) in '", nm, "'."))
  out
}

# ---- Apply explicit type coercions ----
factor_vars <- c(
  "video_file","camera_file","species","site","culture","experiment","well",
  "toxin_exposure","tile","season","collection_site","year",
  "collection_date"
)
for (nm in intersect(factor_vars, names(df))) df[[nm]] <- coerce_factor(df[[nm]], nm)

int_vars <- c("frame_count","mobile_cell_count","stationary_cell_count","total_cells","tile_count","days_from_start")
for (nm in intersect(int_vars, names(df))) df[[nm]] <- coerce_integer(df[[nm]], nm)

num_vars <- c("cu_ug_l","ntu","lux_exposure")
for (nm in intersect(num_vars, names(df))) df[[nm]] <- coerce_numeric(df[[nm]], nm)

date_vars <- c("count_date_yyyy_mm_dd","start_date_yyyy_mm_dd","preparation","release")
for (nm in intersect(date_vars, names(df))) df[[nm]] <- coerce_date_dmy(df[[nm]], nm)

if ("comment" %in% names(df)) {
  if (is.factor(df$comment)) df$comment <- as.character(df$comment)
  df$comment <- as.character(df$comment)
}

# ---- Derived variables (no row drops; log checks) ----
if (all(c("total_cells", "mobile_cell_count") %in% names(df))) {
  df <- df %>% mutate(non_mobile_count = total_cells - mobile_cell_count)
  log_msg("Derived 'non_mobile_count' = total_cells - mobile_cell_count.")
} else {
  log_msg("WARNING: Could not derive non_mobile_count (missing total_cells and/or mobile_cell_count).")
}

if (all(c("total_cells", "mobile_cell_count") %in% names(df))) {
  computed <- ifelse(is.na(df$total_cells) | df$total_cells == 0,
                     NA_real_,
                     df$mobile_cell_count / df$total_cells)
  
  if ("motility_ratio" %in% names(df)) {
    df$motility_ratio <- coerce_numeric(df$motility_ratio, "motility_ratio")
    diff_idx <- which(!is.na(df$motility_ratio) & !is.na(computed) & abs(df$motility_ratio - computed) > 1e-6)
    if (length(diff_idx) > 0) {
      log_msg(paste0(
        "WARNING: motility_ratio differs from mobile/total for ",
        length(diff_idx),
        " row(s). Replacing motility_ratio with computed mobile/total (no rows dropped)."
      ))
    } else {
      log_msg("motility_ratio matches computed mobile/total within tolerance (where comparable).")
    }
  } else {
    log_msg("motility_ratio column missing; creating motility_ratio from mobile/total.")
  }
  
  df$motility_ratio <- computed
  rng <- range(df$motility_ratio, na.rm = TRUE)
  log_msg(paste0("motility_ratio range (na.rm=TRUE): [", signif(rng[1], 4), ", ", signif(rng[2], 4), "]"))
  
  outside <- sum(df$motility_ratio < 0 | df$motility_ratio > 1, na.rm = TRUE)
  if (outside > 0) log_msg(paste0("WARNING: ", outside, " motility_ratio value(s) outside [0,1]. No clamping applied."))
} else {
  log_msg("WARNING: Could not derive motility_ratio (missing total_cells and/or mobile_cell_count).")
}

# ---- Write processed data to parquet (safe replace on Windows) ----
tmp_parq <- paste0(out_parq, ".tmp")

# write to tmp first
arrow::write_parquet(df, tmp_parq)
log_msg(paste0("Wrote temp parquet: ", tmp_parq))

# replace old parquet
if (file.exists(out_parq)) file.remove(out_parq)
file.rename(tmp_parq, out_parq)
log_msg(paste0("Replaced parquet: ", out_parq))


# ---- Write data dictionary ----
dict_spec <- tibble::tribble(
  ~column_name, ~meaning_role, ~unit, ~type_declared, ~notes,
  "video_file", "Unique observation ID", "", "factor", "",
  "camera_file", "File provenance", "", "factor", "",
  "frame_count", "QC check", "", "integer", "",
  "mobile_cell_count", "Raw count – biological response variable", "", "integer", "",
  "stationary_cell_count", "Raw count – biological response variable", "", "integer", "",
  "total_cells", "Raw count – biological response variable", "", "integer", "",
  "motility_ratio", "Derived = mobile/total", "", "numeric (0-1)", "Derived variable; use alongside counts",
  "non_mobile_count", "Derived = total - mobile", "", "integer", "Derived variable",
  "count_date_yyyy_mm_dd", "Date", "", "Date", "dd/mm/yyyy parsed",
  "start_date_yyyy_mm_dd", "Date", "", "Date", "dd/mm/yyyy parsed",
  "species", "Fixed effect", "", "factor", "",
  "site", "Fixed effect", "", "factor", "2-letter site code",
  "culture", "Random effect grouping", "", "factor", "",
  "experiment", "Random effect grouping", "", "factor", "",
  "well", "Experimental unit", "", "factor", "",
  "days_from_start", "Time predictor", "", "integer", "",
  "toxin_exposure", "Treatment indicator", "", "factor", "",
  "cu_ug_l", "Copper concentration", "ug/L", "numeric", "",
  "ntu", "Turbidity", "NTU", "numeric", "",
  "tile", "Surface type", "", "factor", "",
  "tile_count", "Response variable (settlement proxy)", "", "integer", "",
  "lux_exposure", "Light exposure", "lux", "numeric", "",
  "comment", "QC notes", "", "character", "Free text",
  "collection_date", "Metadata", "", "factor", "",
  "year", "Metadata", "", "factor", "",
  "season", "Fixed effect", "", "factor", "Spring/Summer/Autumn/Winter",
  "collection_site", "Fixed effect metadata", "", "factor", "Full site name",
  "preparation", "Metadata date", "", "Date", "",
  "release", "Metadata date", "", "Date", ""
)

actual_types <- tibble(
  column_name = names(df),
  type_actual = vapply(df, function(x) paste(class(x), collapse = "/"), character(1))
)

dict_out <- dict_spec %>%
  left_join(actual_types, by = "column_name") %>%
  mutate(.order = match(column_name, dict_spec$column_name)) %>%
  arrange(.order) %>%
  select(-.order)

readr::write_csv(dict_out, out_dict)
log_msg(paste0("Wrote data dictionary: ", out_dict))

# ---- Toxin exposure codebook + counts (metadata; no data changes) ----
toxin_codebook <- tibble::tribble(
  ~toxin_exposure,      ~toxin_label,                                   ~material_class, ~source,                                            ~notes,
  "CONTROL",            "Control (no added particulate/toxin)",          "control",       NA_character_,                                      NA_character_,
  "none",               "None / unlabelled control",                     "control",       NA_character_,                                      "Kept as-is; consider harmonising to CONTROL later only if you decide.",
  
  "BaSO4",              "Barium sulfate (BaSO₄)",                        "particulate",   "Sigma-Aldrich (243353)",                           "Reference particulate",
  
  "BDC",                "Brake dust (fine; Cambridge HEPA-derived)",     "particulate",   "University of Cambridge (Siriel Saladin)",         "Fine brake wear; hypothesised higher toxicity",
  "BDS",                "Brake dust (coarse; Sussex tyre/lorries)",      "particulate",   "University of Sussex (Phil Howe)",                 "Coarse brake wear; mixed particle sizes",
  
  "Cu0",                "Copper (elemental; Cu⁰)",                       "chemical",      NA_character_,                                      "Speciation test: metallic copper",
  "CuO",                "Copper oxide (CuO)",                            "chemical",      NA_character_,                                      "Speciation test: particulate/oxide form",
  "CuSO4",              "Copper sulfate (CuSO₄; ionic copper)",          "chemical",      NA_character_,                                      "Speciation test: ionic copper",
  
  "SPM",                "Suspended particulate matter (SPM)",            "particulate",   NA_character_,                                      "Field-relevant mixed particulates",
  "SAND",               "Sand",                                          "particulate",   NA_character_,                                      NA_character_,
  "PEAT",               "Peat",                                          "particulate",   NA_character_,                                      NA_character_,
  "GRAPHITE",           "Graphite",                                      "particulate",   NA_character_,                                      NA_character_,
  "PLASTER OF PARIS",   "Plaster of Paris",                              "particulate",   NA_character_,                                      NA_character_,
  "COIR",               "Coir fibres",                                   "particulate",   NA_character_,                                      NA_character_,
  "KAOLINITE",          "Kaolinite clay",                                "particulate",   NA_character_,                                      NA_character_,
  "CLAY SCHOOL",        "Clay (school source)",                          "particulate",   NA_character_,                                      NA_character_,
  "CLAY CLOGGED",       "Clay (clogged / high solids)",                  "particulate",   NA_character_,                                      "As recorded in raw metadata"
)

if ("toxin_exposure" %in% names(df) && is.factor(df$toxin_exposure)) {
  present_codes <- levels(df$toxin_exposure)
  
  toxin_codebook_out <- toxin_codebook %>%
    dplyr::filter(toxin_exposure %in% present_codes)
  
  missing_in_codebook <- setdiff(present_codes, toxin_codebook_out$toxin_exposure)
  if (length(missing_in_codebook) > 0) {
    log_msg(paste0(
      "WARNING: toxin_exposure level(s) not in codebook (left unlabeled): ",
      paste(missing_in_codebook, collapse = ", ")
    ))
  }
  
  out_codebook <- here::here("outputs", "toxin_exposure_codebook.csv")
  readr::write_csv(toxin_codebook_out, out_codebook)
  log_msg(paste0("Wrote toxin exposure codebook: ", out_codebook))
  
  out_counts <- here::here("outputs", "toxin_exposure_counts.csv")
  df %>%
    count(toxin_exposure, name = "n_rows") %>%
    arrange(desc(n_rows)) %>%
    readr::write_csv(out_counts)
  log_msg(paste0("Wrote toxin exposure counts: ", out_counts))
} else {
  log_msg("NOTE: toxin_exposure missing or not a factor; skipped toxin codebook + counts outputs.")
}

# ---- Write data summary ----
col_summary <- tibble(
  column    = names(df),
  class     = vapply(df, function(x) paste(class(x), collapse = "/"), character(1)),
  n_missing = vapply(df, function(x) sum(is.na(x)), integer(1)),
  n_unique  = vapply(df, function(x) dplyr::n_distinct(x, na.rm = FALSE), integer(1))
)

num_preview <- function(x) {
  if (!is.numeric(x)) return(NA_character_)
  if (all(is.na(x))) return(NA_character_)
  r <- range(x, na.rm = TRUE)
  paste0("[", signif(r[1], 5), ", ", signif(r[2], 5), "]")
}

fct_preview <- function(x) {
  if (!is.factor(x)) return(NA_character_)
  lv <- levels(x)
  if (length(lv) <= 8) paste(lv, collapse = ", ")
  else paste0(paste(lv[1:8], collapse = ", "), ", ... (", length(lv), " levels)")
}

col_summary <- col_summary %>%
  mutate(
    range_numeric = vapply(df, num_preview, character(1)),
    levels_factor = vapply(df, fct_preview, character(1))
  )

if (all(c("site","species") %in% names(df))) {
  ss <- df %>% count(site, species, name = "n_rows") %>% arrange(desc(n_rows))
  ss_path <- here::here("outputs","design_counts_site_x_species.csv")
  readr::write_csv(ss, ss_path)
  log_msg(paste0("Wrote design counts: ", ss_path))
}
if ("season" %in% names(df)) {
  se <- df %>% count(season, name = "n_rows") %>% arrange(desc(n_rows))
  se_path <- here::here("outputs","design_counts_season.csv")
  readr::write_csv(se, se_path)
  log_msg(paste0("Wrote design counts: ", se_path))
}

lines <- c(
  "Manuscript 4 — 01_import_clean.R summary",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Input: ", in_path),
  paste0("Rows: ", nrow(df)),
  paste0("Columns: ", ncol(df)),
  "",
  "Notes:",
  "- No rows dropped; no filtering applied.",
  "- motility_ratio recomputed as mobile_cell_count / total_cells (NA if total_cells == 0 or NA).",
  "- non_mobile_count derived as total_cells - mobile_cell_count.",
  "",
  "Column summary:"
)

tab_lines <- capture.output(print(col_summary, n = Inf, width = Inf))
lines <- c(lines, tab_lines, "", "End of summary.")

writeLines(lines, out_sum)
log_msg(paste0("Wrote data summary: ", out_sum))

log_msg("01_import_clean.R complete.")
