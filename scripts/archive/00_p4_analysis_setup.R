# 01_p4_analysis_setup.R
# Purpose: Build an analysis-ready dataset for Paper 4 (SPM + light + brake wear)
# Safe to re-run every session. Produces: data_processed/ms4_p4_analysis.parquet (+ optional CSV)

# ---- 0) Setup ----
source(here::here("scripts", "00_project_setup.R"))

# ---- Packages (explicit; do not rely on global session) ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(stringr)
  library(arrow)
})

log_base("Starting 01_p4_analysis_setup.R")

# ---- 1) Load master cleaned data ----
# Prefer parquet if available (fast + stable types)
master_path_parquet <- here::here("data_processed", "ms4_clean.parquet")
master_path_csv     <- here::here("data_processed", "ms4_clean.csv")  # adjust if yours differs

if (file.exists(master_path_parquet)) {
  log_base(paste("Loading master parquet:", master_path_parquet))
  df <- arrow::read_parquet(master_path_parquet)
} else if (file.exists(master_path_csv)) {
  log_base(paste("Loading master CSV:", master_path_csv))
  df <- readr::read_csv(master_path_csv, show_col_types = FALSE)
} else {
  stop("Could not find master cleaned dataset. Expected either: ",
       master_path_parquet, " or ", master_path_csv)
}

df <- janitor::clean_names(df)

# ---- 1b) Backward compatibility: video_file vs video_file_name ----
# Your canonical column is now video_file (from 01_import_clean.R).
# This block lets older scripts keep working.

# If the parquet has video_file but this script expects video_file_name:
if (!("video_file_name" %in% names(df)) && ("video_file" %in% names(df))) {
  df <- dplyr::mutate(df, video_file_name = as.character(video_file))
}

# If some legacy dataset had video_file_name but not video_file:
if (("video_file_name" %in% names(df)) && !("video_file" %in% names(df))) {
  df <- dplyr::rename(df, video_file = video_file_name)
}

# ---- 2) Minimal sanity checks ----
required_cols <- c(
  "experiment", "species",
  "mobile_cell_count", "stationary_cell_count", "total_cells",
  "toxin_exposure", "ntu", "lux_exposure",
  "well", "video_file"
)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}


# Ensure counts are numeric/integer-ish
df <- df %>%
  mutate(
    mobile_cell_count     = as.numeric(mobile_cell_count),
    stationary_cell_count = as.numeric(stationary_cell_count),
    total_cells           = as.numeric(total_cells),
    ntu                   = suppressWarnings(as.numeric(ntu)),
    lux_exposure           = suppressWarnings(as.numeric(lux_exposure))
  )

# Derive non-mobile count (in case it isn't present)
df <- df %>%
  mutate(
    non_mobile_count = ifelse(!is.na(stationary_cell_count), stationary_cell_count,
                              total_cells - mobile_cell_count),
    motility_ratio = ifelse(total_cells > 0, mobile_cell_count / total_cells, NA_real_)
  )

# Integrity checks
if (any(df$mobile_cell_count > df$total_cells, na.rm = TRUE)) {
  stop("Data integrity error: mobile_cell_count exceeds total_cells in at least one row.")
}
if (any(df$non_mobile_count < 0, na.rm = TRUE)) {
  stop("Data integrity error: non_mobile_count is negative in at least one row.")
}

log_base(paste("Master rows:", nrow(df)))

# ---- 3) Parse experiment block (e.g., 8.2 / 9.2 / 10.2 / 11.2) ----
df <- df %>%
  mutate(
    experiment_chr   = as.character(experiment),
    experiment_block = as.integer(sub("\\..*", "", experiment_chr)),
    experiment_suffix = suppressWarnings(as.integer(sub(".*\\.", "", experiment_chr)))
  )

# ---- 4) Paper 4 subset (experiments 8–11; L. digitata focus) ----
# Adjust this filter if you later include other species.
df_p4 <- df %>%
  filter(experiment_block %in% c(8, 9, 10, 11)) %>%
  filter(species %in% c("L. digitata", "Laminaria digitata", "Ld", "LD") | experiment_suffix == 2) %>%
  mutate(
    species_std = "L. digitata"
  )

log_base(paste("Paper 4 subset rows:", nrow(df_p4)))
log_base(paste("Paper 4 experiments present:", paste(sort(unique(df_p4$experiment_chr)), collapse = ", ")))

# ---- 5) Paper-facing treatment recodes (keep raw codes too) ----
df_p4 <- df_p4 %>%
  mutate(
    toxin_exposure = as.character(toxin_exposure),
    
    # SPM composition (Exp 10.2 mostly)
    spm_material = case_when(
      toxin_exposure %in% c("CONTROL") ~ "Control",
      toxin_exposure %in% c("SPM") ~ "Field SPM",
      toxin_exposure %in% c("PEAT", "COIR") ~ "Organic",
      toxin_exposure %in% c("KAOLINITE", "SAND", "PLASTER OF PARIS", "GRAPHITE",
                            "CLAY CLOGGED", "CLAY SCHOOL") ~ "Mineral / inert",
      TRUE ~ NA_character_
    ),
    
    # Brake wear (Exp 11.2)
    brake_type = case_when(
      toxin_exposure %in% c("none", "CONTROL") ~ "Control",
      toxin_exposure %in% c("BDS") ~ "Brake dust (coarse)",
      toxin_exposure %in% c("BDC") ~ "Brake dust (fine/superfine)",
      toxin_exposure %in% c("BaSO4", "BASO4") ~ "Brake wear (BaSO4)",
      TRUE ~ NA_character_
    ),
    
    brake_class = case_when(
      toxin_exposure %in% c("BDS") ~ "coarse",
      toxin_exposure %in% c("BDC") ~ "fine",
      TRUE ~ NA_character_
    ),
    
    # Convenience transformed predictors
    log_ntu = log10(pmax(ntu, 0) + 1),
    lux_f = factor(lux_exposure, levels = sort(unique(lux_exposure)))
  )

# ---- 6) Create analysis-ready “blocks” as a list (optional) ----
blocks <- list(
  exp8_spm_gradient = df_p4 %>% filter(experiment_chr == "8.2"),
  exp9_light        = df_p4 %>% filter(experiment_chr == "9.2"),
  exp10_composition = df_p4 %>% filter(experiment_chr == "10.2"),
  exp11_brake       = df_p4 %>% filter(experiment_chr == "11.2")
)

# Quick block summaries for the log
for (nm in names(blocks)) {
  d <- blocks[[nm]]
  log_base(paste0(
    nm, ": n=", nrow(d),
    " | toxin_exposure=", paste(sort(unique(d$toxin_exposure)), collapse = ", "),
    " | ntu=", paste(sort(unique(d$ntu)), collapse = ", "),
    " | lux=", paste(sort(unique(d$lux_exposure)), collapse = ", ")
  ))
}

# ---- 7) Save analysis-ready dataset ----
out_parquet <- here::here("data_processed", "ms4_p4_analysis.parquet")
arrow::write_parquet(df_p4, out_parquet)
log_base(paste("Saved:", out_parquet))

# Optional: also save CSV for quick viewing
out_csv <- here::here("data_processed", "ms4_p4_analysis.csv")
readr::write_csv(df_p4, out_csv)
log_base(paste("Saved:", out_csv))

log_base("Finished 00_p4_analysis_setup.R")

