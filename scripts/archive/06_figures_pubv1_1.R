# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      06_figures_pub.R
# Author:      Marianne Glascott
# Date:        2026-02-10
#
# Purpose:
#   Generate publication-ready figures for Manuscript 4 from saved model objects/tables.
#
# Inputs:
#   models/ms4_motility_core_glmmTMB.rds
#   models/ms4_motility_ntu_L_digitata_glmmTMB.rds (optional)
#   data_processed/ms4_clean.parquet (for reference only)
#
# Outputs:
#   figures_pub/*.pdf and figures_pub/*.png
#   outputs/06_figures_pub_log.txt
# ================================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(arrow)
  library(stringr)
  library(ggplot2)
  library(glmmTMB)
  library(emmeans)
  library(rlang)   # for .data pronoun in ggplot
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
dir.create(here::here("figures_pub"), recursive = TRUE, showWarnings = FALSE)

core_rds  <- here::here("models", "ms4_motility_core_glmmTMB.rds")
ntu_rds   <- here::here("models", "ms4_motility_ntu_L_digitata_glmmTMB.rds")
data_parq <- here::here("data_processed", "ms4_clean.parquet")

if (!file.exists(core_rds)) stop("Missing core model: ", core_rds, "\nRun scripts/04_models_motility.R first.")

core_fit <- readRDS(core_rds)
log_msg(paste0("Loaded core model: ", core_rds))

df <- NULL
if (file.exists(data_parq)) {
  df <- arrow::read_parquet(data_parq)
  log_msg(paste0("Loaded parquet for reference: n=", nrow(df), " rows; p=", ncol(df), " columns."))
}

# ---- Plot style helpers ----
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "plain"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "plain"),
      legend.position = "right"
    )
}

save_pub <- function(p, filename_stub, width, height) {
  out_pdf <- here::here("figures_pub", paste0(filename_stub, ".pdf"))
  out_png <- here::here("figures_pub", paste0(filename_stub, ".png"))
  
  # Use base pdf device for maximum Windows compatibility
  ggsave(out_pdf, p, width = width, height = height, units = "in", dpi = 300, device = "pdf")
  ggsave(out_png, p, width = width, height = height, units = "in", dpi = 300)
  
  log_msg(paste0("Saved: ", out_pdf))
  log_msg(paste0("Saved: ", out_png))
}

# ---- Helper: robustly extract response-scale EMMs with CIs ----
emm_to_df <- function(emm_obj) {
  s <- as.data.frame(summary(emm_obj, type = "response", infer = c(TRUE, TRUE)))
  
  # Pick estimate column
  est_col <- intersect(c("response", "prob", "rate", "emmean"), names(s))[1]
  if (is.na(est_col)) stop("Could not find an estimate column in emmeans summary output.")
  
  # Pick CI columns (different naming across models/df assumptions)
  lcl_col <- intersect(c("lower.CL", "asymp.LCL", "lower.HPD"), names(s))[1]
  ucl_col <- intersect(c("upper.CL", "asymp.UCL", "upper.HPD"), names(s))[1]
  if (is.na(lcl_col) || is.na(ucl_col)) {
    stop("Could not find CI columns in emmeans summary output. Columns were: ",
         paste(names(s), collapse = ", "))
  }
  
  s %>%
    rename(
      estimate = !!est_col,
      lcl      = !!lcl_col,
      ucl      = !!ucl_col
    )
}

# ====================================================================
# FIGURE 1: Species × Site (EMMs from core model)
# ====================================================================
emm_ss <- emmeans::emmeans(core_fit, ~ species | site, type = "response")
emm_ss_df <- emm_to_df(emm_ss)

p1 <- ggplot(emm_ss_df, aes(x = species, y = estimate, ymin = lcl, ymax = ucl)) +
  geom_pointrange() +
  facet_wrap(~ site, nrow = 1) +
  labs(
    title = "Zoospore motility by species and site",
    x = NULL,
    y = "Estimated motility (EMM; proportion motile)"
  ) +
  theme_pub()

save_pub(p1, "Fig1_motility_emm_species_by_site", width = 7.0, height = 3.2)

# ====================================================================
# FIGURE 2: Season (EMMs from core model)
# ====================================================================
emm_season <- emmeans::emmeans(core_fit, ~ season, type = "response")
emm_season_df <- emm_to_df(emm_season)

p2 <- ggplot(emm_season_df, aes(x = season, y = estimate, ymin = lcl, ymax = ucl)) +
  geom_pointrange() +
  labs(
    title = "Seasonal variation in zoospore motility",
    x = NULL,
    y = "Estimated motility (EMM; proportion motile)"
  ) +
  theme_pub()

save_pub(p2, "Fig2_motility_emm_season", width = 3.6, height = 3.2)

# ====================================================================
# FIGURE 3: NTU effect in L. digitata (if model exists)
# ====================================================================
if (file.exists(ntu_rds)) {
  
  ntu_fit <- readRDS(ntu_rds)
  log_msg(paste0("Loaded NTU model: ", ntu_rds))
  
  # Use model.frame for correct factor levels (important: experiment_id isn't in parquet)
  mf <- model.frame(ntu_fit)
  
  # NTU range: prefer observed in model frame
  ntu_rng <- range(mf$ntu, na.rm = TRUE)
  ntu_grid <- data.frame(ntu = seq(ntu_rng[1], ntu_rng[2], length.out = 120))
  
  # Reference values for other covariates (from model frame for consistency)
  site_ref   <- levels(mf$site)[1]
  season_ref <- levels(mf$season)[1]
  
  tox_levels <- levels(mf$toxin_exposure)
  tox_ref <- if ("CONTROL" %in% tox_levels) "CONTROL" else tox_levels[1]
  
  lux_ref <- if ("lux_exposure" %in% names(mf)) median(mf$lux_exposure, na.rm = TRUE) else 0
  cu_ref  <- if ("cu_ug_l" %in% names(mf)) median(mf$cu_ug_l, na.rm = TRUE) else 0
  
  culture_ref <- levels(mf$culture)[1]
  expid_ref   <- levels(mf$experiment_id)[1]
  
  newdata <- ntu_grid %>%
    mutate(
      site           = factor(site_ref,   levels = levels(mf$site)),
      season         = factor(season_ref, levels = levels(mf$season)),
      toxin_exposure = factor(tox_ref,    levels = levels(mf$toxin_exposure)),
      lux_exposure   = lux_ref,
      cu_ug_l        = cu_ref,
      culture        = factor(culture_ref, levels = levels(mf$culture)),
      experiment_id  = factor(expid_ref,   levels = levels(mf$experiment_id))
    )
  
  # Predictions on LINK scale + transform (stable CI)
  pr <- predict(
    ntu_fit,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    allow.new.levels = TRUE
  )
  
  pred_df <- newdata %>%
    transmute(
      ntu = ntu,
      fit = plogis(pr$fit),
      lwr = plogis(pr$fit - 1.96 * pr$se.fit),
      upr = plogis(pr$fit + 1.96 * pr$se.fit)
    )
  
  p3 <- ggplot(pred_df, aes(x = ntu, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line() +
    labs(
      title = "Effect of turbidity (NTU) on motility in L. digitata",
      x = "NTU",
      y = "Predicted motility (proportion motile)"
    ) +
    theme_pub()
  
  save_pub(p3, "Fig3_motility_ntu_L_digitata", width = 4.6, height = 3.2)
  
} else {
  log_msg("NOTE: NTU model RDS not found; skipping Fig3. (Expected only if identifiable.)")
}

log_msg("06_figures_pub.R complete.")
