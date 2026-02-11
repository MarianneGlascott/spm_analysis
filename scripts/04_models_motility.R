# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      04_models_motility.R
# Author:      Marianne Glascott
# Date:        2026-02-10
#
# Purpose:
#   Fit pre-specified mixed-effects models for Manuscript 4 (motility).
#   Save ALL model objects, diagnostics, and tables to disk.
#
# Inputs:
#   data_processed/ms4_clean.parquet
#
# Outputs (required, written to disk):
#   models/ms4_motility_core_glmmTMB.rds
#   models/ms4_motility_ntu_<species>_glmmTMB.rds (0+ files; only if identifiable)
#   tables/model_core_fixed_effects.csv
#   tables/model_core_emm_species_by_site.csv
#   tables/model_core_emm_season.csv
#   tables/model_ntu_fixed_effects_<species>.csv
#   outputs/model_summary_core.txt
#   outputs/model_summary_ntu_<species>.txt
#   outputs/model_diagnostics_core.txt
#   outputs/model_diagnostics_ntu_<species>.txt
#   outputs/04_models_motility_log.txt
#   outputs/experiment_id_map.csv
#   outputs/design_site_x_species_counts.csv
#   outputs/model_rows_dropped_core.csv
#   outputs/model_rows_dropped_ntu_<species>.csv
#
# Rules:
#   - No silent filtering: any NA-based row drops are logged and written to disk
#   - No invented QC thresholds
#   - Respect design: site × species crossed; season fixed
#   - Random effects: (1 | culture / experiment_id)
#   - Primary response: motility_ratio
# ================================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(arrow)
  library(stringr)
  library(tidyr)
  
  # Modelling
  library(glmmTMB)
  library(emmeans)
  library(broom.mixed)
  library(DHARMa)
})

# ---- Logging ----
log_file <- here::here("outputs", "04_models_motility_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 04_models_motility.R")

# ---- Paths ----
in_path <- here::here("data_processed", "ms4_clean.parquet")
dir.create(here::here("models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_path)) {
  stop("Input parquet not found at: ", in_path, "\nRun scripts/01_import_clean.R first.")
}

df <- arrow::read_parquet(in_path)
log_msg(paste0("Loaded parquet: n=", nrow(df), " rows; p=", ncol(df), " columns."))

# ---- Guardrail checks (no dropping) ----
required <- c("motility_ratio","species","site","season","culture","experiment",
              "ntu","lux_exposure","cu_ug_l","toxin_exposure")
missing_req <- setdiff(required, names(df))
if (length(missing_req) > 0) {
  stop("Missing required column(s) for modelling: ", paste(missing_req, collapse = ", "))
}

# Use video_file as row id if present; else create one
if (!("video_file" %in% names(df))) {
  df <- df %>% mutate(video_file = paste0("row_", row_number()))
  log_msg("NOTE: video_file not present; created synthetic row id 'row_<n>'.")
}

# ---- Option A: experiment_id for missing experiments (culture-specific controls) ----
df <- df %>%
  mutate(
    experiment_chr = as.character(experiment),
    culture_chr    = as.character(culture),
    experiment_id  = if_else(is.na(experiment_chr) | experiment_chr == "",
                             paste0("CTRL_", culture_chr),
                             experiment_chr)
  )

# Write audit map (unique mappings)
exp_map <- df %>%
  distinct(culture, experiment, experiment_id) %>%
  arrange(culture, experiment_id)

out_map <- here::here("outputs", "experiment_id_map.csv")
readr::write_csv(exp_map, out_map)
log_msg(paste0("Wrote experiment_id map: ", out_map))

n_ctrl <- sum(is.na(df$experiment) | df$experiment == "")
log_msg(paste0("Rows with experiment missing: ", n_ctrl, " (encoded as CTRL_<culture>)."))
log_msg(paste0("Distinct experiment_id levels: ", dplyr::n_distinct(df$experiment_id)))

# Ensure factors
df <- df %>%
  mutate(
    species        = factor(species),
    site           = factor(site),
    season         = factor(season),
    culture        = factor(culture),
    experiment_id  = factor(experiment_id),
    toxin_exposure = factor(toxin_exposure)
  )

# ---- Design completeness audit (site × species) ----
design_ss <- df %>%
  count(site, species, name = "n_rows") %>%
  tidyr::complete(site, species, fill = list(n_rows = 0)) %>%
  arrange(site, species)

design_path <- here::here("outputs", "design_site_x_species_counts.csv")
readr::write_csv(design_ss, design_path)
log_msg(paste0("Wrote design completeness table: ", design_path))

missing_cells <- design_ss %>% filter(n_rows == 0)
if (nrow(missing_cells) > 0) {
  log_msg(paste0(
    "Design is incomplete; missing site×species cells: ",
    paste(paste0(missing_cells$site, "×", missing_cells$species), collapse = "; "),
    ". This can cause rank-deficient interaction terms in species*site."
  ))
}

# ---- Beta squeeze for motility_ratio ----
squeeze_beta <- function(y) {
  n <- sum(!is.na(y))
  if (n == 0) return(y)
  (y * (n - 1) + 0.5) / n
}
df <- df %>% mutate(motility_beta = squeeze_beta(motility_ratio))

n0 <- sum(df$motility_ratio == 0, na.rm = TRUE)
n1 <- sum(df$motility_ratio == 1, na.rm = TRUE)
log_msg(paste0("motility_ratio exact zeros: ", n0, "; exact ones: ", n1, " (handled via beta squeeze)."))

# ---- Stabilising transforms (do NOT overwrite raw columns) ----
# copper spans huge range: use log10(x+1) + scale to improve optimiser stability
df <- df %>%
  mutate(
    cu_log10 = log10(cu_ug_l + 1),
    cu_z     = as.numeric(scale(cu_log10)),
    lux_z    = as.numeric(scale(lux_exposure)),
    ntu_z    = as.numeric(scale(ntu))
  )

log_msg("Created transformed/scaled covariates: cu_log10, cu_z, lux_z, ntu_z (raw columns retained).")

# ---- Helper: prepare model data WITHOUT silent NA dropping ----
prepare_model_data <- function(formula, data, model_name) {
  vars <- unique(all.vars(formula))
  vars <- setdiff(vars, c("motility_beta")) # present in data but included anyway
  needed <- intersect(vars, names(data))
  
  df_sub <- data[, c("video_file", needed), drop = FALSE]
  cc <- stats::complete.cases(df_sub)
  
  dropped <- data %>%
    mutate(.complete = cc) %>%
    filter(!.complete) %>%
    select(video_file, all_of(needed))
  
  used <- data %>% filter(cc)
  
  out_drop <- here::here("outputs", paste0("model_rows_dropped_", model_name, ".csv"))
  readr::write_csv(dropped, out_drop)
  
  log_msg(paste0(
    "Model '", model_name, "': complete cases = ", nrow(used),
    " / ", nrow(data),
    " (dropped ", nrow(dropped), " rows with NA in model terms). Wrote: ", out_drop
  ))
  
  list(data_used = used, dropped = dropped, dropped_path = out_drop)
}

# ---- Helper: fit model and capture warnings into log ----
fit_glmmTMB_logged <- function(formula, data, model_name) {
  warn <- character(0)
  
  w_handler <- function(w) {
    warn <<- c(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  
  fit <- withCallingHandlers(
    glmmTMB::glmmTMB(
      formula = formula,
      family  = beta_family(link = "logit"),
      data    = data,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 1e4, eval.max = 1e4)
      )
    ),
    warning = w_handler
  )
  
  if (length(warn) > 0) {
    log_msg(paste0("Model '", model_name, "' warnings (", length(warn), "):"))
    for (i in seq_along(warn)) log_msg(paste0("  [", i, "] ", warn[i]))
  } else {
    log_msg(paste0("Model '", model_name, "': no warnings captured."))
  }
  
  list(fit = fit, warnings = warn)
}

# ================================================================================
# CORE MODEL (all data)
# ================================================================================
core_formula <- motility_beta ~
  species * site + season +
  lux_z + cu_z + toxin_exposure +
  (1 | culture / experiment_id)

log_msg("Fitting CORE model (beta/logit): species*site + season + lux_z + cu_z + toxin + (1|culture/experiment_id)")

core_prep <- prepare_model_data(core_formula, df, model_name = "core")
core_res  <- fit_glmmTMB_logged(core_formula, core_prep$data_used, model_name = "core")
core_fit  <- core_res$fit

# ---- CORE: convergence check (logged + saved) ----
core_pdHess <- isTRUE(core_fit$sdr$pdHess)
core_conv   <- core_fit$fit$convergence

log_msg(paste0("CORE pdHess (TRUE is good): ", core_pdHess))
log_msg(paste0("CORE optimizer convergence code (0 is good): ", core_conv))

if (!core_pdHess) {
  log_msg("WARNING: Core model pdHess is FALSE (non-positive-definite Hessian). Sensitivity models in Script 05 will help identify a stable specification.")
}

core_sum_path <- here::here("outputs", "model_summary_core.txt")
writeLines(capture.output(summary(core_fit)), core_sum_path)
log_msg(paste0("Wrote core model summary: ", core_sum_path))

core_rds <- here::here("models", "ms4_motility_core_glmmTMB.rds")
saveRDS(core_fit, core_rds)
log_msg(paste0("Saved core model: ", core_rds))

# ---- CORE: Fixed effects table ----
core_fix <- broom.mixed::tidy(core_fit, effects = "fixed", conf.int = TRUE) %>%
  arrange(term)

out_core_fix <- here::here("tables", "model_core_fixed_effects.csv")
readr::write_csv(core_fix, out_core_fix)
log_msg(paste0("Wrote core fixed effects table: ", out_core_fix))

# ---- CORE: EMMs ----
emm_ss <- emmeans::emmeans(core_fit, ~ species | site, type = "response")
emm_ss_df <- as.data.frame(emm_ss)

out_emm_ss <- here::here("tables", "model_core_emm_species_by_site.csv")
readr::write_csv(emm_ss_df, out_emm_ss)
log_msg(paste0("Wrote core EMMs species|site: ", out_emm_ss))

emm_season <- emmeans::emmeans(core_fit, ~ season, type = "response")
emm_season_df <- as.data.frame(emm_season)

out_emm_season <- here::here("tables", "model_core_emm_season.csv")
readr::write_csv(emm_season_df, out_emm_season)
log_msg(paste0("Wrote core EMMs season: ", out_emm_season))

# ---- CORE diagnostics (DHARMa) ----
diag_core_path <- here::here("outputs", "model_diagnostics_core.txt")
sim_core <- DHARMa::simulateResiduals(core_fit, n = 500)

diag_lines <- c(
  "Core model diagnostics (DHARMa)",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Model:",
  deparse(core_formula),
  "",
  paste0("pdHess: ", core_pdHess),
  paste0("optimizer convergence code: ", core_conv),
  "",
  "DHARMa tests:"
)
diag_lines <- c(diag_lines, capture.output(print(DHARMa::testDispersion(sim_core))))
diag_lines <- c(diag_lines, capture.output(print(DHARMa::testZeroInflation(sim_core))))
diag_lines <- c(diag_lines, capture.output(print(DHARMa::testUniformity(sim_core))))

writeLines(diag_lines, diag_core_path)
log_msg(paste0("Wrote core model diagnostics: ", diag_core_path))

# ================================================================================
# NTU MODELS (species-specific; only if identifiable)
# ================================================================================
species_levels <- levels(df$species)

for (sp in species_levels) {
  
  df_sp <- df %>% filter(species == sp)
  
  ntu_unique <- dplyr::n_distinct(df_sp$ntu[!is.na(df_sp$ntu)])
  log_msg(paste0("Species ", sp, ": n=", nrow(df_sp), "; unique NTU (non-missing) = ", ntu_unique))
  
  if (ntu_unique < 2) {
    log_msg(paste0("Skipping NTU slope model for ", sp, " (NTU not identifiable: <2 unique values)."))
    next
  }
  
  sp_slug <- sp %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_$", "") %>%
    str_replace_all("^_", "")
  
  ntu_formula <- motility_beta ~
    ntu_z + site + season + lux_z + cu_z + toxin_exposure +
    (1 | culture / experiment_id)
  
  log_msg(paste0("Fitting NTU model for ", sp, " (beta/logit): ntu_z + site + season + lux_z + cu_z + toxin + (1|culture/experiment_id)"))
  
  prep_sp <- prepare_model_data(ntu_formula, df_sp, model_name = paste0("ntu_", sp_slug))
  res_sp  <- fit_glmmTMB_logged(ntu_formula, prep_sp$data_used, model_name = paste0("ntu_", sp_slug))
  fit_sp  <- res_sp$fit
  
  sp_pdHess <- isTRUE(fit_sp$sdr$pdHess)
  sp_conv   <- fit_sp$fit$convergence
  
  log_msg(paste0("NTU pdHess for ", sp, " (TRUE is good): ", sp_pdHess))
  log_msg(paste0("NTU optimizer convergence code for ", sp, " (0 is good): ", sp_conv))
  
  if (!sp_pdHess) {
    log_msg(paste0("WARNING: NTU model pdHess FALSE for ", sp, ". Script 05 will test alternative specs if needed."))
  }
  
  sp_sum_path <- here::here("outputs", paste0("model_summary_ntu_", sp_slug, ".txt"))
  writeLines(capture.output(summary(fit_sp)), sp_sum_path)
  log_msg(paste0("Wrote NTU model summary: ", sp_sum_path))
  
  rds_path <- here::here("models", paste0("ms4_motility_ntu_", sp_slug, "_glmmTMB.rds"))
  saveRDS(fit_sp, rds_path)
  log_msg(paste0("Saved NTU model: ", rds_path))
  
  fix_sp <- broom.mixed::tidy(fit_sp, effects = "fixed", conf.int = TRUE) %>%
    arrange(term)
  
  out_fix_sp <- here::here("tables", paste0("model_ntu_fixed_effects_", sp_slug, ".csv"))
  readr::write_csv(fix_sp, out_fix_sp)
  log_msg(paste0("Wrote NTU fixed effects table: ", out_fix_sp))
  
  diag_path <- here::here("outputs", paste0("model_diagnostics_ntu_", sp_slug, ".txt"))
  sim_sp <- DHARMa::simulateResiduals(fit_sp, n = 500)
  
  lines_sp <- c(
    paste0("NTU model diagnostics (", sp, ")"),
    paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Model:",
    deparse(ntu_formula),
    "",
    paste0("pdHess: ", sp_pdHess),
    paste0("optimizer convergence code: ", sp_conv),
    "",
    "DHARMa tests:"
  )
  
  lines_sp <- c(lines_sp, capture.output(print(DHARMa::testDispersion(sim_sp))))
  lines_sp <- c(lines_sp, capture.output(print(DHARMa::testZeroInflation(sim_sp))))
  lines_sp <- c(lines_sp, capture.output(print(DHARMa::testUniformity(sim_sp))))
  
  writeLines(lines_sp, diag_path)
  log_msg(paste0("Wrote NTU model diagnostics: ", diag_path))
}

log_msg("04_models_motility.R complete.")
