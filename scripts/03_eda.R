# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      03_eda.R
# Author:      Marianne Glascott
# Date:        2026-02-10
#
# Purpose:
#   Exploratory data analysis (EDA) for Manuscript 4.
#   Descriptive diagnostics only — NO modelling, NO filtering.
#
# Input:
#   data_processed/ms4_clean.parquet
#
# Outputs (written to disk):
#   figures/eda_*.png
#   tables/eda_*.csv
#   outputs/eda_summary.txt
#   outputs/03_eda_log.txt
#
# Rules:
#   - No row removal
#   - No transformations chosen here
#   - No statistical tests or models
#   - Motility_ratio is the primary endpoint
# ================================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(arrow)
  library(tidyr)
  library(forcats)
})

# ---- Logging ----
log_file <- here::here("outputs", "03_eda_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 03_eda.R")

# ---- Paths ----
in_path <- here::here("data_processed", "ms4_clean.parquet")

fig_dir <- here::here("figures")
tab_dir <- here::here("tables")
out_sum <- here::here("outputs", "eda_summary.txt")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_path)) {
  stop("Input parquet not found at: ", in_path, "\nRun scripts/01_import_clean.R first.")
}

df <- arrow::read_parquet(in_path)
log_msg(paste0("Loaded parquet: n=", nrow(df), " rows; p=", ncol(df), " columns."))

# ---- EDA 1: Design balance tables ----

design_counts <- df %>%
  count(site, species, season, name = "n_rows") %>%
  arrange(site, species, season)

readr::write_csv(
  design_counts,
  here::here("tables", "eda_design_counts_site_species_season.csv")
)
log_msg("Wrote design balance table: site × species × season.")

# ---- EDA 2: Distribution of motility_ratio (overall) ----

p_overall <- ggplot(df, aes(x = motility_ratio)) +
  geom_histogram(bins = 40, colour = "grey30", fill = "grey80") +
  labs(
    title = "Distribution of zoospore motility ratio (all observations)",
    x = "Motility ratio (mobile / total)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("figures", "eda_motility_ratio_overall.png"),
  plot = p_overall,
  width = 7,
  height = 5,
  dpi = 300
)
log_msg("Saved overall motility_ratio distribution.")

# ---- EDA 3: Motility by species ----

p_species <- ggplot(df, aes(x = species, y = motility_ratio)) +
  geom_violin(fill = "grey85", colour = "grey40", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    title = "Zoospore motility ratio by species",
    x = "Species",
    y = "Motility ratio"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("figures", "eda_motility_ratio_by_species.png"),
  plot = p_species,
  width = 7,
  height = 5,
  dpi = 300
)
log_msg("Saved motility_ratio by species.")

# ---- EDA 4: Motility by site × species ----

p_site_species <- ggplot(df, aes(x = site, y = motility_ratio)) +
  geom_violin(fill = "grey85", colour = "grey40", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ species) +
  labs(
    title = "Zoospore motility ratio by site and species",
    x = "Site",
    y = "Motility ratio"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("figures", "eda_motility_ratio_by_site_species.png"),
  plot = p_site_species,
  width = 9,
  height = 5,
  dpi = 300
)
log_msg("Saved motility_ratio by site × species.")

# ---- EDA 5: Seasonal patterns ----

p_season <- ggplot(df, aes(x = season, y = motility_ratio)) +
  geom_violin(fill = "grey85", colour = "grey40", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    title = "Zoospore motility ratio by season",
    x = "Season",
    y = "Motility ratio"
  ) +
  theme_minimal()

ggsave(
  filename = here::here("figures", "eda_motility_ratio_by_season.png"),
  plot = p_season,
  width = 7,
  height = 5,
  dpi = 300
)
log_msg("Saved motility_ratio by season.")

# ---- EDA 6: Motility vs environmental gradients (scatter only) ----

env_vars <- c("ntu", "lux_exposure", "cu_ug_l")

for (v in env_vars) {
  if (v %in% names(df)) {
    
    df_plot <- df %>%
      filter(!is.na(.data[[v]]), !is.na(motility_ratio))
    
    n_dropped <- nrow(df) - nrow(df_plot)
    log_msg(paste0("Plotting motility_ratio vs ", v, ": omitted ", n_dropped,
                   " row(s) due to NA in ", v, " and/or motility_ratio (plotting only)."))
    
    p <- ggplot(df_plot, aes(x = .data[[v]], y = motility_ratio)) +
      geom_point(alpha = 0.4, size = 1) +
      facet_wrap(~ species) +
      labs(
        title = paste("Motility ratio vs", v),
        x = v,
        y = "Motility ratio"
      ) +
      theme_minimal()
    
    ggsave(
      filename = here::here("figures", paste0("eda_motility_ratio_vs_", v, ".png")),
      plot = p,
      width = 9,
      height = 5,
      dpi = 300
    )
    log_msg(paste0("Saved motility_ratio vs ", v, "."))
  }
}


# ---- EDA 7: Missingness snapshot (descriptive) ----

missing_summary <- df %>%
  summarise(across(
    everything(),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  arrange(desc(n_missing))

readr::write_csv(
  missing_summary,
  here::here("tables", "eda_missingness_summary.csv")
)
log_msg("Wrote missingness summary table.")

# ---- Write EDA summary text ----

lines <- c(
  "Manuscript 4 — 03_eda.R summary",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Rows: ", nrow(df)),
  paste0("Columns: ", ncol(df)),
  "",
  "EDA scope:",
  "- Descriptive only",
  "- No modelling",
  "- No filtering or exclusions",
  "- Motility_ratio treated as primary endpoint",
  "",
  "Key outputs:",
  "- Distribution of motility_ratio (overall, by species, site × species, season)",
  "- Scatterplots of motility_ratio vs NTU, lux_exposure, cu_ug_l",
  "- Design balance table (site × species × season)",
  "- Missingness summary"
)

writeLines(lines, out_sum)
log_msg("Wrote EDA summary text.")

log_msg("03_eda.R complete.")

