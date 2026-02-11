# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      07_figures_refine.R
# Author:      Marianne Glascott
# Date:        2026-02-11
#
# Purpose:
#   Iterative, step-by-step figure development (publication quality).
#   Run ONE figure at a time while you refine.
#
# Inputs:
#   models/ms4_motility_core_glmmTMB.rds
#   models/ms4_motility_ntu_L_digitata_glmmTMB.rds  (optional)
#   data_processed/ms4_clean.parquet
#
# Outputs:
#   figures_pub/*.pdf + *.png (only when save = TRUE)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(arrow)
  library(ggplot2)
  library(stringr)
  library(emmeans)
  library(glmmTMB)
  library(forcats)
})

# ---- Paths ----
fig_dir <- here::here("figures_pub")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

core_path <- here::here("models", "ms4_motility_core_glmmTMB.rds")
ntu_ld_path <- here::here("models", "ms4_motility_ntu_L_digitata_glmmTMB.rds")
dat_path  <- here::here("data_processed", "ms4_clean.parquet")

stopifnot(file.exists(core_path), file.exists(dat_path))

m_core <- readRDS(core_path)
df     <- arrow::read_parquet(dat_path)

# ---- kelp species palette ----
# IMPORTANT: keep names EXACTLY as in df$species levels
species_pal <- c(
  "L. digitata"    = "#c09c0e",
  "L. hyperborea"  = "#687fa2",
  "S. latissima"   = "#ba89c3"
)

# ---- Minimal publication theme (feel free to swap with your existing theme) ----
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "right",
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.minor = element_blank()
    )
}

# ---- Save helper (pdf + png, consistent sizes) ----
save_fig <- function(p, filename_stub, w = 170, h = 120, dpi = 300) {
  # w/h in mm. 170mm ~ 1.5 columns depending on journal
  pdf_path <- file.path(fig_dir, paste0(filename_stub, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename_stub, ".png"))
  
  ggsave(pdf_path, plot = p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(png_path, plot = p, width = w, height = h, units = "mm", dpi = dpi)
  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

# ==============================================================================
# FIG 1 — EMM species by site (core model)
# ==============================================================================

fig1_species_by_site <- function(save = FALSE) {
  emm <- emmeans::emmeans(m_core, ~ species | site, type = "response", infer = c(TRUE, TRUE))
  emm_df <- as.data.frame(emm)
  
  # You currently have asymp.LCL / asymp.UCL (not lower.CL)
  p <- ggplot(emm_df, aes(x = species, y = response, colour = species)) +
    geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                    linewidth = 0.5, fatten = 2) +
    facet_wrap(~ site, nrow = 1) +
    scale_colour_manual(values = species_pal, drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                       limits = c(0, 1)) +
    labs(x = NULL, y = "Motility ratio (EMM, response scale)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (save) save_fig(p, "Fig1_motility_emm_species_by_site", w = 180, h = 110)
  p
}

# ==============================================================================
# FIG 2 — EMM season (core model)  [yours was blank — this guards against that]
# ==============================================================================

fig2_season <- function(save = FALSE) {
  emm <- emmeans::emmeans(m_core, ~ season, type = "response", infer = c(TRUE, TRUE))
  emm_df <- as.data.frame(emm)
  
  # Ensure season is ordered (edit order if you prefer)
  emm_df <- emm_df %>%
    mutate(season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")))
  
  # If this ever goes blank, it's usually because limits drop all points.
  # So: no hard limits unless we confirm range.
  yr <- range(c(emm_df$asymp.LCL, emm_df$asymp.UCL), na.rm = TRUE)
  pad <- 0.03
  ylims <- c(max(0, yr[1] - pad), min(1, yr[2] + pad))
  
  p <- ggplot(emm_df, aes(x = season, y = response)) +
    geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                    linewidth = 0.5, colour = "black", fatten = 2) +
    scale_y_continuous(limits = ylims,
                       labels = scales::label_number(accuracy = 0.01)) +
    labs(x = NULL, y = "Motility ratio (EMM, response scale)") +
    theme_pub()
  
  if (save) save_fig(p, "Fig2_motility_emm_season", w = 140, h = 110)
  p
}

# ==============================================================================
# Helper: build newdata grid for a fitted model (robust to z-scored covariates)
# ==============================================================================

# This works even if your NTU model uses ntu_z / lux_z / cu_z etc.
# We hold other covariates at 0 (their mean) when z-scored, or at median if raw.
make_newdata_for_term <- function(m, term, grid_n = 50) {
  mf <- model.frame(m)
  
  # Identify covariates present in model frame
  vars <- names(mf)
  
  # Create baseline newdata with one row per site/season/toxin if present
  # (keep it simple first; you can expand later)
  newdata <- tibble::tibble()
  
  if ("site" %in% vars)   newdata$site <- levels(mf$site)[1]
  if ("season" %in% vars) newdata$season <- levels(mf$season)[1]
  if ("toxin_exposure" %in% vars) newdata$toxin_exposure <- levels(mf$toxin_exposure)[1]
  
  # Random effects grouping
  if ("culture" %in% vars) newdata$culture <- levels(mf$culture)[1]
  if ("experiment_id" %in% vars) newdata$experiment_id <- levels(mf$experiment_id)[1]
  
  # For all numeric covariates in mf, set baseline at 0 if looks z-scored, else median
  num_vars <- vars[vapply(mf, is.numeric, logical(1))]
  for (v in setdiff(num_vars, c("motility_beta", term))) {
    if (str_detect(v, "_z$")) {
      newdata[[v]] <- 0
    } else {
      newdata[[v]] <- median(mf[[v]], na.rm = TRUE)
    }
  }
  
  # Now create the focal sequence for `term`
  x <- mf[[term]]
  if (!is.numeric(x)) stop("Requested term is not numeric in model frame: ", term)
  
  seq_x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = grid_n)
  newdata <- newdata[rep(1, grid_n), , drop = FALSE]
  newdata[[term]] <- seq_x
  
  newdata
}

# Predict on response scale for glmmTMB beta(logit)
predict_response <- function(m, newdata) {
  eta <- predict(m, newdata = newdata, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  tibble::tibble(
    x = newdata[[setdiff(names(newdata), names(newdata))[1]]]
  )
}

# ==============================================================================
# FIG 3 — L. digitata NTU effect (from NTU model)
# ==============================================================================

fig3_ld_ntu <- function(save = FALSE, show_points = TRUE) {
  if (!file.exists(ntu_ld_path)) stop("NTU model not found: ", ntu_ld_path)
  m_ntu <- readRDS(ntu_ld_path)
  mf <- model.frame(m_ntu)
  
  # Determine term name used in model (ntu_z vs ntu)
  term <- if ("ntu_z" %in% names(mf)) "ntu_z" else if ("ntu" %in% names(mf)) "ntu" else stop("No ntu term found.")
  
  nd <- make_newdata_for_term(m_ntu, term = term, grid_n = 80)
  
  # Predict mean + CI on response scale
  pr <- predict(m_ntu, newdata = nd, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  nd$fit_link <- pr$fit
  nd$se_link  <- pr$se.fit
  
  invlogit <- function(z) 1/(1+exp(-z))
  nd <- nd %>%
    mutate(
      fit = invlogit(fit_link),
      lcl = invlogit(fit_link - 1.96*se_link),
      ucl = invlogit(fit_link + 1.96*se_link)
    )
  
  # For display axis: if z-scored, we can show z OR convert back.
  # For now: label as z if *_z.
  xlab <- if (term == "ntu_z") "Turbidity (NTU, z-scored)" else "Turbidity (NTU)"
  
  p <- ggplot(nd, aes(x = .data[[term]], y = fit)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
    geom_line(linewidth = 0.8) +
    labs(x = xlab, y = "Predicted motility ratio") +
    theme_pub()
  
  if (show_points) {
    # Overlay raw points (from parquet) for context
    df_ld <- df %>% filter(species == "L. digitata") %>% filter(!is.na(ntu), !is.na(motility_ratio))
    # If model is z-scored, map x to z-score used by model frame if available
    if (term == "ntu_z") {
      mu <- mean(df_ld$ntu, na.rm = TRUE)
      sdv <- sd(df_ld$ntu, na.rm = TRUE)
      df_ld <- df_ld %>% mutate(ntu_z = (ntu - mu)/sdv)
      p <- p + geom_point(data = df_ld, aes(x = ntu_z, y = motility_ratio), alpha = 0.25)
    } else {
      p <- p + geom_point(data = df_ld, aes(x = ntu, y = motility_ratio), alpha = 0.25)
    }
  }
  
  if (save) save_fig(p, "Fig3_Ldigitata_motility_vs_ntu", w = 160, h = 120)
  p
}

# ==============================================================================
# FIG 4 — L. digitata LUX effect (from same NTU model)
# ==============================================================================

fig4_ld_lux <- function(save = FALSE, show_points = TRUE) {
  if (!file.exists(ntu_ld_path)) stop("NTU model not found: ", ntu_ld_path)
  m_ntu <- readRDS(ntu_ld_path)
  mf <- model.frame(m_ntu)
  
  term <- if ("lux_z" %in% names(mf)) "lux_z" else if ("lux_exposure" %in% names(mf)) "lux_exposure" else stop("No lux term found.")
  nd <- make_newdata_for_term(m_ntu, term = term, grid_n = 80)
  
  pr <- predict(m_ntu, newdata = nd, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  nd$fit_link <- pr$fit
  nd$se_link  <- pr$se.fit
  
  invlogit <- function(z) 1/(1+exp(-z))
  nd <- nd %>%
    mutate(
      fit = invlogit(fit_link),
      lcl = invlogit(fit_link - 1.96*se_link),
      ucl = invlogit(fit_link + 1.96*se_link)
    )
  
  xlab <- if (term == "lux_z") "Light (lux, z-scored)" else "Light (lux)"
  
  p <- ggplot(nd, aes(x = .data[[term]], y = fit)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
    geom_line(linewidth = 0.8) +
    labs(x = xlab, y = "Predicted motility ratio") +
    theme_pub()
  
  if (show_points) {
    df_ld <- df %>% filter(species == "L. digitata") %>% filter(!is.na(lux_exposure), !is.na(motility_ratio))
    if (term == "lux_z") {
      mu <- mean(df_ld$lux_exposure, na.rm = TRUE)
      sdv <- sd(df_ld$lux_exposure, na.rm = TRUE)
      df_ld <- df_ld %>% mutate(lux_z = (lux_exposure - mu)/sdv)
      p <- p + geom_point(data = df_ld, aes(x = lux_z, y = motility_ratio), alpha = 0.25)
    } else {
      p <- p + geom_point(data = df_ld, aes(x = lux_exposure, y = motility_ratio), alpha = 0.25)
    }
  }
  
  if (save) save_fig(p, "Fig4_Ldigitata_motility_vs_lux", w = 160, h = 120)
  p
}

# ==============================================================================
# DEV RUN SECTION — run lines ONE AT A TIME in the console
# ==============================================================================

# fig1_species_by_site(save = FALSE)
# fig2_season(save = FALSE)
# fig3_ld_ntu(save = FALSE)
# fig4_ld_lux(save = FALSE)

# Once you're happy:
# fig1_species_by_site(save = TRUE)
# fig2_season(save = TRUE)
# fig3_ld_ntu(save = TRUE)
# fig4_ld_lux(save = TRUE)

