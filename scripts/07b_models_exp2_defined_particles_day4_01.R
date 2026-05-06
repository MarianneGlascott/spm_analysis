# ==============================================================================
# Script title: 07b_models_exp2_defined_particles_day4_01.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Fit and visualise Day 4 models for Experiment 2,
#          testing whether defined particle treatments alter
#          kelp zoospore motility across nominal turbidity
#          levels.
# Inputs: Derived Day 4 analysis dataset:
#         /data_derived/day4/ms4_day4_analysis_derived.rds
#         Helper scripts in /R where available
# Outputs: Experiment 2 Day 4 analysis subset; model objects;
#          model comparison tables; fixed-effect summaries;
#          estimated marginal means; treatment contrasts;
#          diagnostic outputs; quick-view and model-estimated
#          figures; plain-text log and result note.
# Date created: 06 May 2026
# Last updated: 06 May 2026
# Notes/dependencies:
# - This script builds on the Day 4 derived dataset created
#   by the preceding data-derivation workflow.
# - The script uses strict expected column names from the
#   derived Day 4 dataset and fails early if the data contract
#   changes.
# - Rows are not dropped except where explicit exclusion flags
#   have already been assigned upstream, and where model-specific
#   complete-case requirements are applied.
# - The primary response is motile fraction, modelled from
#   motile and non-motile cell counts where available.
# - The main biological question is whether motility differs
#   among defined particle types and/or nominal NTU treatments
#   on Day 4.
# - The Type-III ANOVA table is optional and is only produced
#   if the car package is installed; estimated marginal means
#   and treatment contrasts should remain the main interpretive
#   outputs.
# ==============================================================================
# -----------------------------------------------------------------------------
# 0. Housekeeping
# -----------------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(glmmTMB)
  library(DHARMa)
  library(performance)
  library(broom.mixed)
  library(emmeans)
  library(car)
  library(patchwork)
  library(scales)
  library(readr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Paths and helper functions
# -----------------------------------------------------------------------------

# Expected project root: C:/SPM_Analysis
# If running from the RStudio project, here::here() should resolve correctly.

paths <- list(
  data_derived = here("data_derived"),
  outputs      = here("outputs"),
  models       = here("models"),
  tables       = here("tables"),
  figures      = here("figures"),
  logs         = here("outputs", "logs")
)

walk(paths, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# Load helper functions if present. These are optional, so the script can still run
# if one helper file has not yet been created.
helper_files <- c(
  here("R", "helpers_theme.R"),
  here("R", "helpers_save_figures.R"),
  here("R", "helpers_labels.R"),
  here("R", "helpers_model_checks.R"),
  here("R", "helpers_tables.R")
)

existing_helpers <- helper_files[file.exists(helper_files)]
walk(existing_helpers, source)

# Fallback plotting theme if helpers_theme.R has not defined theme_ms4().
if (!exists("theme_ms4")) {
  theme_ms4 <- function(base_size = 11) {
    theme_classic(base_size = base_size) +
      theme(
        axis.title = element_text(size = base_size),
        axis.text  = element_text(size = base_size - 1),
        plot.title = element_text(face = "bold", size = base_size + 1),
        plot.subtitle = element_text(size = base_size),
        legend.title = element_text(size = base_size - 1),
        legend.text  = element_text(size = base_size - 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )
  }
}

# Fallback figure saver if helpers_save_figures.R has not defined save_ms4_figure().
if (!exists("save_ms4_figure")) {
  save_ms4_figure <- function(plot, filename_stub, width = 180, height = 120, units = "mm", dpi = 600) {
    ggsave(
      filename = file.path(paths$figures, paste0(filename_stub, ".pdf")),
      plot = plot,
      width = width,
      height = height,
      units = units,
      device = cairo_pdf
    )
    ggsave(
      filename = file.path(paths$figures, paste0(filename_stub, ".png")),
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi
    )
  }
}

script_id <- "07b_models_exp2_defined_particles_day4_01"
log_file <- file.path(paths$logs, paste0(script_id, "_log.txt"))

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

writeLines(c(
  paste0("Script: ", script_id),
  paste0("Run time: ", Sys.time()),
  paste0("Working directory: ", getwd()),
  ""
), con = log_file)

# -----------------------------------------------------------------------------
# 2. Load derived Day 4 dataset
# -----------------------------------------------------------------------------

# Primary expected file from the previous Day 4 workflow.
derived_rds <- file.path(paths$data_derived, "ms4_day4_analysis_derived.rds")
derived_csv <- file.path(paths$data_derived, "ms4_day4_analysis_derived.csv")

if (file.exists(derived_rds)) {
  day4_all <- readRDS(derived_rds)
  log_msg("Loaded Day 4 dataset from: ", derived_rds)
} else if (file.exists(derived_csv)) {
  day4_all <- read_csv(derived_csv, show_col_types = FALSE)
  log_msg("Loaded Day 4 dataset from: ", derived_csv)
} else {
  stop("Could not find derived Day 4 dataset at either:\n",
       derived_rds, "\n", derived_csv)
}

log_msg("Rows: ", nrow(day4_all))
log_msg("Columns: ", ncol(day4_all))

# -----------------------------------------------------------------------------
# 3. Column checks and flexible variable mapping
# -----------------------------------------------------------------------------

# The script allows for small differences in naming between earlier scripts.
# Edit this section only if your current dataset uses different names.

candidate_cols <- list(
  experiment       = c("experiment", "experiment_number", "exp"),
  species          = c("species"),
  site             = c("site", "collection_site"),
  culture          = c("culture", "culture_id"),
  well             = c("well", "well_id"),
  treatment        = c("particle_treatment", "particle_type", "particle", "toxin_exposure", "treatment"),
  ntu              = c("ntu", "ntu_exposure", "target_ntu", "nominal_ntu"),
  motility_ratio   = c("motility_ratio", "motile_fraction", "mobile_fraction", "prop_motile"),
  mobile_count     = c("mobile_cell_count", "motile_count", "mobile_count"),
  stationary_count = c("stationary_cell_count", "non_mobile_count", "immobile_count", "stationary_count"),
  total_count      = c("total_cells", "total_cell_count", "n_cells", "cell_count"),
  days_from_start  = c("days_from_start", "day", "culture_day")
)

pick_col <- function(data, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Missing required column for ", label %||% paste(candidates, collapse = "/"),
         ". Tried: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

`%||%` <- function(x, y) if (is.null(x)) y else x

col_experiment       <- pick_col(day4_all, candidate_cols$experiment, required = FALSE, label = "experiment")
col_species          <- pick_col(day4_all, candidate_cols$species, required = FALSE, label = "species")
col_site             <- pick_col(day4_all, candidate_cols$site, required = FALSE, label = "site")
col_culture          <- pick_col(day4_all, candidate_cols$culture, required = FALSE, label = "culture")
col_well             <- pick_col(day4_all, candidate_cols$well, required = TRUE,  label = "well")
col_treatment        <- pick_col(day4_all, candidate_cols$treatment, required = TRUE, label = "particle treatment")
col_ntu              <- pick_col(day4_all, candidate_cols$ntu, required = TRUE, label = "NTU")
col_motility_ratio   <- pick_col(day4_all, candidate_cols$motility_ratio, required = FALSE, label = "motility ratio")
col_mobile_count     <- pick_col(day4_all, candidate_cols$mobile_count, required = FALSE, label = "mobile count")
col_stationary_count <- pick_col(day4_all, candidate_cols$stationary_count, required = FALSE, label = "stationary/non-mobile count")
col_total_count      <- pick_col(day4_all, candidate_cols$total_count, required = FALSE, label = "total count")
col_days_from_start  <- pick_col(day4_all, candidate_cols$days_from_start, required = FALSE, label = "days from start")

log_msg("\nColumn mapping:")
log_msg("  experiment       = ", col_experiment)
log_msg("  species          = ", col_species)
log_msg("  site             = ", col_site)
log_msg("  culture          = ", col_culture)
log_msg("  well             = ", col_well)
log_msg("  treatment        = ", col_treatment)
log_msg("  ntu              = ", col_ntu)
log_msg("  motility_ratio   = ", col_motility_ratio)
log_msg("  mobile_count     = ", col_mobile_count)
log_msg("  stationary_count = ", col_stationary_count)
log_msg("  total_count      = ", col_total_count)
log_msg("  days_from_start  = ", col_days_from_start)

# -----------------------------------------------------------------------------
# 4. Create Experiment 2 Day 4 analysis subset
# -----------------------------------------------------------------------------

exp2_day4 <- day4_all

# Filter to Experiment 2 if an experiment column exists.
if (!is.na(col_experiment)) {
  exp2_day4 <- exp2_day4 %>%
    filter(.data[[col_experiment]] %in% c(2, "2", "Exp2", "Experiment 2", "defined_particles"))
}

# Confirm or derive Day 4 if needed. Because this file should already be Day 4,
# this only acts as a safeguard where a day column is present.
if (!is.na(col_days_from_start)) {
  possible_day_values <- unique(exp2_day4[[col_days_from_start]])
  if (any(possible_day_values == 4, na.rm = TRUE)) {
    exp2_day4 <- exp2_day4 %>% filter(.data[[col_days_from_start]] == 4)
  }
}

# Standardised analysis columns.
exp2_day4 <- exp2_day4 %>%
  mutate(
    particle_type = factor(.data[[col_treatment]]),
    ntu_numeric   = as.numeric(.data[[col_ntu]]),
    ntu_factor    = factor(ntu_numeric, levels = sort(unique(ntu_numeric))),
    well_id       = factor(.data[[col_well]])
  )

if (!is.na(col_species)) {
  exp2_day4 <- exp2_day4 %>% mutate(species = factor(.data[[col_species]]))
} else {
  exp2_day4 <- exp2_day4 %>% mutate(species = factor("not_recorded"))
}

if (!is.na(col_site)) {
  exp2_day4 <- exp2_day4 %>% mutate(site = factor(.data[[col_site]]))
} else {
  exp2_day4 <- exp2_day4 %>% mutate(site = factor("not_recorded"))
}

if (!is.na(col_culture)) {
  exp2_day4 <- exp2_day4 %>% mutate(culture_id = factor(.data[[col_culture]]))
} else {
  exp2_day4 <- exp2_day4 %>% mutate(culture_id = factor("not_recorded"))
}

# Derive motile/non-motile counts where possible.
if (!is.na(col_mobile_count) && !is.na(col_stationary_count)) {
  exp2_day4 <- exp2_day4 %>%
    mutate(
      motile_n     = as.integer(.data[[col_mobile_count]]),
      non_motile_n = as.integer(.data[[col_stationary_count]]),
      total_n      = motile_n + non_motile_n,
      motile_prop  = motile_n / total_n
    )
} else if (!is.na(col_mobile_count) && !is.na(col_total_count)) {
  exp2_day4 <- exp2_day4 %>%
    mutate(
      motile_n     = as.integer(.data[[col_mobile_count]]),
      total_n      = as.integer(.data[[col_total_count]]),
      non_motile_n = total_n - motile_n,
      motile_prop  = motile_n / total_n
    )
} else if (!is.na(col_motility_ratio) && !is.na(col_total_count)) {
  exp2_day4 <- exp2_day4 %>%
    mutate(
      total_n      = as.integer(.data[[col_total_count]]),
      motile_prop  = as.numeric(.data[[col_motility_ratio]]),
      motile_n     = round(motile_prop * total_n),
      non_motile_n = total_n - motile_n
    )
} else if (!is.na(col_motility_ratio)) {
  exp2_day4 <- exp2_day4 %>%
    mutate(
      motile_prop = as.numeric(.data[[col_motility_ratio]]),
      total_n = NA_integer_,
      motile_n = NA_integer_,
      non_motile_n = NA_integer_
    )
} else {
  stop("Could not derive motility response. Need either counts or motility ratio.")
}

# Basic checks.
if (nrow(exp2_day4) == 0) {
  stop("Experiment 2 Day 4 subset has zero rows. Check experiment/day labels.")
}

summary_tbl <- exp2_day4 %>%
  summarise(
    rows = n(),
    wells = n_distinct(well_id),
    species_n = n_distinct(species),
    cultures_n = n_distinct(culture_id),
    particle_types = paste(sort(unique(as.character(particle_type))), collapse = ", "),
    ntu_levels = paste(sort(unique(ntu_numeric)), collapse = ", "),
    min_motile_prop = min(motile_prop, na.rm = TRUE),
    max_motile_prop = max(motile_prop, na.rm = TRUE)
  )

write_csv(summary_tbl, file.path(paths$tables, paste0(script_id, "_dataset_summary.csv")))

log_msg("\nExperiment 2 Day 4 subset created.")
log_msg("Rows: ", nrow(exp2_day4))
log_msg("Wells: ", n_distinct(exp2_day4$well_id))
log_msg("Particle types: ", paste(sort(unique(as.character(exp2_day4$particle_type))), collapse = ", "))
log_msg("NTU levels: ", paste(sort(unique(exp2_day4$ntu_numeric)), collapse = ", "))
log_msg("Cultures: ", n_distinct(exp2_day4$culture_id))

# Save analysis subset for transparency.
saveRDS(exp2_day4, file.path(paths$data_derived, paste0(script_id, "_analysis_subset.rds")))
write_csv(exp2_day4, file.path(paths$data_derived, paste0(script_id, "_analysis_subset.csv")))

# -----------------------------------------------------------------------------
# 5. Exploratory summaries and quick-view plot
# -----------------------------------------------------------------------------

summary_by_treatment <- exp2_day4 %>%
  group_by(particle_type, ntu_factor) %>%
  summarise(
    n_wells = n(),
    mean_motile_prop = mean(motile_prop, na.rm = TRUE),
    sd_motile_prop = sd(motile_prop, na.rm = TRUE),
    se_motile_prop = sd_motile_prop / sqrt(n_wells),
    median_motile_prop = median(motile_prop, na.rm = TRUE),
    min_motile_prop = min(motile_prop, na.rm = TRUE),
    max_motile_prop = max(motile_prop, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_by_treatment, file.path(paths$tables, paste0(script_id, "_summary_by_particle_ntu.csv")))

p_raw <- ggplot(exp2_day4, aes(x = ntu_factor, y = motile_prop)) +
  geom_point(aes(shape = particle_type),
             position = position_jitter(width = 0.08, height = 0),
             alpha = 0.8,
             size = 2) +
  stat_summary(aes(group = particle_type),
               fun = mean,
               geom = "line",
               linewidth = 0.6) +
  stat_summary(aes(group = particle_type),
               fun = mean,
               geom = "point",
               size = 2.6) +
  facet_wrap(~ particle_type, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Experiment 2 Day 4: motile fraction by defined particle treatment",
    subtitle = "Points are wells; lines show treatment means across NTU levels",
    x = "Nominal turbidity treatment (NTU)",
    y = "Motile fraction",
    shape = "Particle type"
  ) +
  theme_ms4()

save_ms4_figure(p_raw, paste0(script_id, "_raw_quick_view"), width = 210, height = 115)

# -----------------------------------------------------------------------------
# 6. Model-ready data
# -----------------------------------------------------------------------------

# For binomial models, retain rows with complete counts.
model_dat_counts <- exp2_day4 %>%
  filter(
    !is.na(motile_n),
    !is.na(non_motile_n),
    !is.na(total_n),
    total_n > 0,
    motile_n >= 0,
    non_motile_n >= 0,
    !is.na(particle_type),
    !is.na(ntu_factor)
  ) %>%
  mutate(
    particle_type = droplevels(particle_type),
    ntu_factor = droplevels(ntu_factor),
    species = droplevels(species),
    site = droplevels(site),
    culture_id = droplevels(culture_id)
  )

# For beta/proportion fallback models, retain valid proportions away from exact 0/1.
model_dat_prop <- exp2_day4 %>%
  filter(
    !is.na(motile_prop),
    motile_prop >= 0,
    motile_prop <= 1,
    !is.na(particle_type),
    !is.na(ntu_factor)
  ) %>%
  mutate(
    motile_prop_beta = (motile_prop * (n() - 1) + 0.5) / n(),
    particle_type = droplevels(particle_type),
    ntu_factor = droplevels(ntu_factor),
    species = droplevels(species),
    site = droplevels(site),
    culture_id = droplevels(culture_id)
  )

log_msg("\nModel-ready rows with counts: ", nrow(model_dat_counts))
log_msg("Model-ready rows with proportions: ", nrow(model_dat_prop))

# -----------------------------------------------------------------------------
# 7. Fit candidate models
# -----------------------------------------------------------------------------

models <- list()
model_notes <- tibble(model = character(), note = character())

safe_fit <- function(model_name, expr) {
  log_msg("\nFitting model: ", model_name)
  fit <- tryCatch(
    expr,
    error = function(e) {
      log_msg("  FAILED: ", conditionMessage(e))
      model_notes <<- bind_rows(model_notes, tibble(model = model_name, note = conditionMessage(e)))
      NULL
    },
    warning = function(w) {
      log_msg("  WARNING: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  fit
}

# Main model: binomial GLMM where counts are available.
# Random intercept for culture is included only if there is more than one culture.
use_culture_re <- n_distinct(model_dat_counts$culture_id) > 1

if (nrow(model_dat_counts) > 0) {
  if (use_culture_re) {
    models$binom_interaction <- safe_fit(
      "binom_interaction",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type * ntu_factor + (1 | culture_id),
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_additive <- safe_fit(
      "binom_additive",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type + ntu_factor + (1 | culture_id),
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_ntu_only <- safe_fit(
      "binom_ntu_only",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ ntu_factor + (1 | culture_id),
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_particle_only <- safe_fit(
      "binom_particle_only",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type + (1 | culture_id),
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
  } else {
    models$binom_interaction <- safe_fit(
      "binom_interaction",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type * ntu_factor,
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_additive <- safe_fit(
      "binom_additive",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type + ntu_factor,
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_ntu_only <- safe_fit(
      "binom_ntu_only",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ ntu_factor,
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
    
    models$binom_particle_only <- safe_fit(
      "binom_particle_only",
      glmmTMB(
        cbind(motile_n, non_motile_n) ~ particle_type,
        family = binomial(link = "logit"),
        data = model_dat_counts
      )
    )
  }
}

# Fallback/sensitivity model: beta regression on well-level proportions.
# This is useful if the count response is unavailable or as a sensitivity check.
if (nrow(model_dat_prop) > 0) {
  use_culture_re_prop <- n_distinct(model_dat_prop$culture_id) > 1
  
  if (use_culture_re_prop) {
    models$beta_interaction <- safe_fit(
      "beta_interaction",
      glmmTMB(
        motile_prop_beta ~ particle_type * ntu_factor + (1 | culture_id),
        family = beta_family(link = "logit"),
        data = model_dat_prop
      )
    )
  } else {
    models$beta_interaction <- safe_fit(
      "beta_interaction",
      glmmTMB(
        motile_prop_beta ~ particle_type * ntu_factor,
        family = beta_family(link = "logit"),
        data = model_dat_prop
      )
    )
  }
}

models <- models[!map_lgl(models, is.null)]

if (length(models) == 0) {
  stop("No models fitted successfully. Check response columns and treatment structure.")
}

# Save model objects.
saveRDS(models, file.path(paths$models, paste0(script_id, "_models.rds")))

# -----------------------------------------------------------------------------
# 8. Model comparison and summaries
# -----------------------------------------------------------------------------

model_comparison <- map_dfr(names(models), function(nm) {
  fit <- models[[nm]]
  tibble(
    model = nm,
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    df = attr(logLik(fit), "df")
  )
}) %>%
  arrange(AIC)

write_csv(model_comparison, file.path(paths$tables, paste0(script_id, "_model_comparison_aic.csv")))
write_csv(model_notes, file.path(paths$tables, paste0(script_id, "_model_notes.csv")))

log_msg("\nModel comparison written to tables.")
log_msg("Best AIC model: ", model_comparison$model[1])

# Prefer binomial model for inference where available; otherwise use beta model.
primary_model_name <- case_when(
  "binom_interaction" %in% names(models) ~ "binom_interaction",
  "binom_additive" %in% names(models) ~ "binom_additive",
  TRUE ~ model_comparison$model[1]
)

primary_model <- models[[primary_model_name]]
log_msg("Primary model selected for estimated means: ", primary_model_name)

# Fixed-effect summary.
fixed_effects <- tidy(primary_model, effects = "fixed", conf.int = TRUE)
write_csv(fixed_effects, file.path(paths$tables, paste0(script_id, "_primary_model_fixed_effects.csv")))

# Type-II/III tests. For interaction models this gives a useful broad check, but
# interpretation should focus on estimated marginal means and contrasts.
anova_primary <- tryCatch({
  car::Anova(primary_model, type = 3) %>%
    as.data.frame() %>%
    rownames_to_column("term")
}, error = function(e) {
  tibble(term = "Anova failed", message = conditionMessage(e))
})

write_csv(anova_primary, file.path(paths$tables, paste0(script_id, "_primary_model_anova_type3.csv")))

# -----------------------------------------------------------------------------
# 9. Estimated marginal means and contrasts
# -----------------------------------------------------------------------------

emm_particle_ntu <- emmeans(primary_model, ~ particle_type * ntu_factor, type = "response")
emm_particle_ntu_tbl <- as_tibble(emm_particle_ntu)
write_csv(emm_particle_ntu_tbl, file.path(paths$tables, paste0(script_id, "_emmeans_particle_by_ntu_response.csv")))

# Within-particle NTU contrasts.
contrasts_ntu_within_particle <- tryCatch({
  contrast(emm_particle_ntu, method = "pairwise", by = "particle_type", adjust = "tukey") %>%
    as_tibble()
}, error = function(e) {
  tibble(note = paste("NTU within particle contrasts failed:", conditionMessage(e)))
})
write_csv(contrasts_ntu_within_particle, file.path(paths$tables, paste0(script_id, "_contrasts_ntu_within_particle.csv")))

# Within-NTU particle contrasts.
contrasts_particle_within_ntu <- tryCatch({
  contrast(emm_particle_ntu, method = "pairwise", by = "ntu_factor", adjust = "tukey") %>%
    as_tibble()
}, error = function(e) {
  tibble(note = paste("Particle within NTU contrasts failed:", conditionMessage(e)))
})
write_csv(contrasts_particle_within_ntu, file.path(paths$tables, paste0(script_id, "_contrasts_particle_within_ntu.csv")))

# -----------------------------------------------------------------------------
# 10. Model diagnostics
# -----------------------------------------------------------------------------

diagnostics_dir <- file.path(paths$outputs, "diagnostics", script_id)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

# DHARMa residual diagnostics.
png(file.path(diagnostics_dir, paste0(primary_model_name, "_dharma_residuals.png")),
    width = 1800, height = 1400, res = 200)
res_primary <- simulateResiduals(primary_model, n = 1000)
plot(res_primary)
dev.off()

# Save DHARMa tests.
dharma_tests <- tibble(
  test = c("uniformity", "dispersion", "outliers"),
  p_value = c(
    tryCatch(testUniformity(res_primary)$p.value, error = function(e) NA_real_),
    tryCatch(testDispersion(res_primary)$p.value, error = function(e) NA_real_),
    tryCatch(testOutliers(res_primary)$p.value, error = function(e) NA_real_)
  )
)
write_csv(dharma_tests, file.path(paths$tables, paste0(script_id, "_primary_model_dharma_tests.csv")))

# performance checks.
perf_check <- tryCatch({
  performance::check_model(primary_model)
}, error = function(e) NULL)

if (!is.null(perf_check)) {
  png(file.path(diagnostics_dir, paste0(primary_model_name, "_performance_check_model.png")),
      width = 1800, height = 1400, res = 200)
  print(perf_check)
  dev.off()
}

# -----------------------------------------------------------------------------
# 11. Prediction plot: observed + model-estimated means
# -----------------------------------------------------------------------------

emm_plot_dat <- emm_particle_ntu_tbl %>%
  rename(
    estimated_motile_prop = prob,
    lower.CL = asymp.LCL,
    upper.CL = asymp.UCL
  )

# The emmeans response column can differ by family; handle common alternatives.
if (!"estimated_motile_prop" %in% names(emm_plot_dat)) {
  response_col <- intersect(c("response", "rate", "emmean"), names(emm_plot_dat))[1]
  if (!is.na(response_col)) {
    emm_plot_dat <- emm_plot_dat %>% rename(estimated_motile_prop = all_of(response_col))
  }
}

if (!"lower.CL" %in% names(emm_plot_dat)) {
  lcl_col <- intersect(c("asymp.LCL", "lower.CL", "LCL"), names(emm_plot_dat))[1]
  if (!is.na(lcl_col)) emm_plot_dat <- emm_plot_dat %>% rename(lower.CL = all_of(lcl_col))
}

if (!"upper.CL" %in% names(emm_plot_dat)) {
  ucl_col <- intersect(c("asymp.UCL", "upper.CL", "UCL"), names(emm_plot_dat))[1]
  if (!is.na(ucl_col)) emm_plot_dat <- emm_plot_dat %>% rename(upper.CL = all_of(ucl_col))
}

p_model <- ggplot() +
  geom_point(
    data = exp2_day4,
    aes(x = ntu_factor, y = motile_prop),
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.45,
    size = 1.8
  ) +
  geom_errorbar(
    data = emm_plot_dat,
    aes(x = ntu_factor, ymin = lower.CL, ymax = upper.CL, group = particle_type),
    width = 0.08,
    linewidth = 0.4
  ) +
  geom_line(
    data = emm_plot_dat,
    aes(x = ntu_factor, y = estimated_motile_prop, group = particle_type),
    linewidth = 0.7
  ) +
  geom_point(
    data = emm_plot_dat,
    aes(x = ntu_factor, y = estimated_motile_prop),
    size = 2.8
  ) +
  facet_wrap(~ particle_type, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Experiment 2 Day 4: defined particles and zoospore motility",
    subtitle = paste0("Observed wells with model-estimated means and 95% confidence intervals; model: ", primary_model_name),
    x = "Nominal turbidity treatment (NTU)",
    y = "Motile fraction"
  ) +
  theme_ms4()

save_ms4_figure(p_model, paste0(script_id, "_model_estimated_means"), width = 210, height = 115)

# -----------------------------------------------------------------------------
# 12. Optional compact manuscript-style result note
# -----------------------------------------------------------------------------

result_note <- c(
  paste0("Experiment 2 Day 4 defined-particle model summary"),
  paste0("Run time: ", Sys.time()),
  "",
  paste0("Analysis subset rows: ", nrow(exp2_day4)),
  paste0("Particle types: ", paste(sort(unique(as.character(exp2_day4$particle_type))), collapse = ", ")),
  paste0("NTU levels: ", paste(sort(unique(exp2_day4$ntu_numeric)), collapse = ", ")),
  paste0("Primary model: ", primary_model_name),
  "",
  "Model comparison by AIC:",
  capture.output(print(model_comparison)),
  "",
  "Primary model fixed effects:",
  capture.output(print(fixed_effects)),
  "",
  "Type-III test table:",
  capture.output(print(anova_primary)),
  "",
  "DHARMa tests:",
  capture.output(print(dharma_tests)),
  "",
  "Interpretation note:",
  "Use the estimated marginal means and treatment contrasts rather than raw coefficient signs for narrative interpretation, especially if the particle_type × NTU interaction is retained. The key biological question is whether motile fraction changes consistently with increasing NTU within each particle type, or whether particle identity alters the response independently of nominal turbidity."
)

writeLines(result_note, file.path(paths$outputs, paste0(script_id, "_result_note.txt")))

# -----------------------------------------------------------------------------
# 13. End log
# -----------------------------------------------------------------------------

log_msg("\nOutputs written:")
log_msg("  Analysis subset: ", file.path(paths$data_derived, paste0(script_id, "_analysis_subset.rds")))
log_msg("  Models: ", file.path(paths$models, paste0(script_id, "_models.rds")))
log_msg("  Tables: ", paths$tables)
log_msg("  Figures: ", paths$figures)
log_msg("  Diagnostics: ", diagnostics_dir)
log_msg("  Result note: ", file.path(paths$outputs, paste0(script_id, "_result_note.txt")))
log_msg("\nScript complete.")

sessionInfo_file <- file.path(paths$logs, paste0(script_id, "_sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), sessionInfo_file)
