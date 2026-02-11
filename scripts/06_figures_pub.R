# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      06_figures_pub.R
# Author:      Marianne Glascott
# Date:        2026-02-10
#
# Purpose:
#   Publication-ready figures for Manuscript 4 (motility)
#
# Inputs:
#   models/ms4_motility_core_glmmTMB.rds
#   models/ms4_motility_ntu_<species>_glmmTMB.rds (optional)
#   data_processed/ms4_clean.parquet
#
# Outputs:
#   figures_pub/*.png
#   outputs/06_figures_pub_log.txt
# ================================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(arrow)
  library(stringr)
  library(ggplot2)
  library(emmeans)
  library(glmmTMB)
})

# ---- Logging ----
log_file <- here::here("outputs", "06_figures_pub_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 06_figures_pub.R")

# ---- Paths ----
model_path <- here::here("models", "ms4_motility_core_glmmTMB.rds")
data_path  <- here::here("data_processed", "ms4_clean.parquet")

fig_dir <- here::here("figures_pub")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_path)) stop("Core model not found at: ", model_path)
if (!file.exists(data_path))  stop("Processed parquet not found at: ", data_path)

core_fit <- readRDS(model_path)
log_msg(paste0("Loaded core model: ", normalizePath(model_path, winslash = "/")))

df <- arrow::read_parquet(data_path)
log_msg(paste0("Loaded parquet for reference: n=", nrow(df), " rows; p=", ncol(df), " columns."))

# ---- Styling helpers ----
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

save_pub <- function(p, filename, width = 8, height = 4.5, dpi = 300) {
  out <- file.path(fig_dir, filename)
  ggsave(out, plot = p, width = width, height = height, dpi = dpi, units = "in")
  log_msg(paste0("Saved figure: ", normalizePath(out, winslash = "/")))
}

# ---- EMM helper: standardise to columns: estimate, lwr, upr ----
emm_to_df <- function(emm_obj) {
  
  s <- as.data.frame(summary(emm_obj, type = "response", infer = c(TRUE, TRUE)))
  
  # Identify estimate column on response scale
  est_col <- if ("response" %in% names(s)) "response" else {
    cand <- intersect(names(s), c("prob", "rate", "emmean"))
    if (length(cand) == 0) stop("Could not identify response column in emmeans summary output. Names: ",
                                paste(names(s), collapse = ", "))
    cand[1]
  }
  
  # Identify CI columns across emmeans versions
  lcl_cand <- intersect(names(s), c("lower.CL", "asymp.LCL", "LCL", "lower.HPD"))
  ucl_cand <- intersect(names(s), c("upper.CL", "asymp.UCL", "UCL", "upper.HPD"))
  
  if (length(lcl_cand) == 0 || length(ucl_cand) == 0) {
    stop("EMMeans summary did not include CI columns. Names were: ",
         paste(names(s), collapse = ", "),
         "\nTry updating emmeans OR ensure infer=c(TRUE,TRUE) works in your version.")
  }
  
  s <- s %>%
    mutate(
      estimate = .data[[est_col]],
      lwr = .data[[lcl_cand[1]]],
      upr = .data[[ucl_cand[1]]]
    )
  
  s
}

#===============================================================================
# FIGURE 1: Species x Site EMMs
#===============================================================================
emm_ss <- emmeans::emmeans(core_fit, ~ species | site, type = "response")
emm_ss_df <- emm_to_df(emm_ss)

p1 <- ggplot(
  emm_ss_df,
  aes(x = species, y = estimate, ymin = lwr, ymax = upr)
) +
  geom_pointrange() +
  facet_wrap(~ site, nrow = 1) +
  labs(
    title = "Zoospore motility by species and site",
    x = NULL,
    y = "Estimated motility (EMM; proportion motile)"
  ) +
  theme_pub()

save_pub(p1, "Fig1_emm_species_by_site.png", width = 10, height = 4)

#===============================================================================
# FIGURE 2: Season EMMs
#===============================================================================
emm_season <- emmeans::emmeans(core_fit, ~ season, type = "response")
emm_se_df <- emm_to_df(emm_season)

p2 <- ggplot(
  emm_se_df,
  aes(x = season, y = estimate, ymin = lwr, ymax = upr)
) +
  geom_pointrange() +
  labs(
    title = "Seasonal variation in zoospore motility",
    x = NULL,
    y = "Estimated motility (EMM; proportion motile)"
  ) +
  theme_pub()

save_pub(p2, "Fig2_emm_season.png", width = 7, height = 4)

#===============================================================================
# FIGURE 3: NTU effect in L. digitata (if model exists)
#===============================================================================
ntu_rds <- here::here("models", "ms4_motility_ntu_L_digitata_glmmTMB.rds")

if (file.exists(ntu_rds)) {
  
  ntu_fit <- readRDS(ntu_rds)
  log_msg(paste0("Loaded NTU model: ", normalizePath(ntu_rds, winslash = "/")))
  
  df_ld <- df %>% filter(species == "L. digitata" & !is.na(ntu))
  ntu_rng <- range(df_ld$ntu, na.rm = TRUE)
  ntu_grid <- data.frame(ntu = seq(ntu_rng[1], ntu_rng[2], length.out = 100))
  
  # Reference values (interpretable curve)
  site_ref   <- levels(factor(df_ld$site))[1]
  season_ref <- levels(factor(df_ld$season))[1]
  tox_levels <- levels(factor(df_ld$toxin_exposure))
  tox_ref    <- if ("CONTROL" %in% tox_levels) "CONTROL" else tox_levels[1]
  
  lux_ref <- median(df_ld$lux_exposure, na.rm = TRUE)
  cu_ref  <- median(df_ld$cu_ug_l, na.rm = TRUE)
  
  # Need valid factor levels for random effects columns used in 04
  culture_ref <- levels(factor(df_ld$culture))[1]
  expid_ref <- {
    expid <- ifelse(is.na(df_ld$experiment) | df_ld$experiment == "",
                    paste0("CTRL_", as.character(df_ld$culture)),
                    as.character(df_ld$experiment))
    levels(factor(expid))[1]
  }
  
  newdata <- ntu_grid %>%
    mutate(
      site = factor(site_ref, levels = levels(factor(df_ld$site))),
      season = factor(season_ref, levels = levels(factor(df_ld$season))),
      toxin_exposure = factor(tox_ref, levels = levels(factor(df_ld$toxin_exposure))),
      lux_exposure = lux_ref,
      cu_ug_l = cu_ref,
      culture = factor(culture_ref, levels = levels(factor(df_ld$culture))),
      experiment_id = factor(expid_ref)
    )
  
  # CI on link scale -> back-transform
  pr_link <- predict(ntu_fit, newdata = newdata, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  
  pred_df <- newdata %>%
    mutate(
      eta = pr_link$fit,
      se  = pr_link$se.fit,
      lwr = plogis(eta - 1.96 * se),
      fit = plogis(eta),
      upr = plogis(eta + 1.96 * se)
    ) %>%
    select(ntu, fit, lwr, upr)
  
  p3 <- ggplot(pred_df, aes(x = ntu, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line() +
    labs(
      title = "Effect of turbidity (NTU) on motility in L. digitata",
      x = "NTU",
      y = "Predicted motility (proportion motile)"
    ) +
    theme_pub()
  
  save_pub(p3, "Fig3_motility_ntu_L_digitata.png", width = 6, height = 4)
  
} else {
  log_msg("NOTE: NTU model RDS not found; skipping Fig3.")
}

log_msg("06_figures_pub.R complete.")



