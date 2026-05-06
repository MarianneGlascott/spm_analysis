# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      05_models_sensitivity_motility.R
# Author:      Marianne Glascott
# Date:        2026-02-10
#
# Purpose:
#   Sensitivity checks for Manuscript 4 motility models to support robustness:
#     - Address rank-deficiency in species×site interaction (incomplete crossing)
#     - Check convergence / pdHess stability
#     - Check dependence on random-effects structure
#     - Provide model comparison table + diagnostics summaries
#
# Inputs:
#   data_processed/ms4_clean.parquet
#
# Outputs (written to disk):
#   models/ms4_motility_sens_*.rds
#   outputs/05_models_sensitivity_motility_log.txt
#   outputs/sensitivity_model_comparison.csv
#   outputs/sensitivity_fixed_effects_all.csv
#   outputs/sensitivity_summaries.txt
#   outputs/model_summary_sens_*.txt
#   outputs/model_diagnostics_sens_*.txt
#   outputs/model_rows_dropped_sens_*.csv
#
# Rules:
#   - No silent filtering: any NA-based row drops are logged and written
#   - No invented QC thresholds
#   - Sensitivity models are complementary; primary inference remains Script 04
# ================================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(arrow)
  library(stringr)
  library(tidyr)
  
  library(glmmTMB)
  library(broom.mixed)
  library(DHARMa)
  library(tibble)
})

# ---- Logging ----
log_file <- here::here("outputs", "05_models_sensitivity_motility_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

log_msg("Starting 05_models_sensitivity_motility.R")

# ---- Paths ----
in_path <- here::here("data_processed", "ms4_clean.parquet")
dir.create(here::here("models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_path)) stop("Input parquet not found at: ", in_path, "\nRun scripts/01_import_clean.R first.")

df <- arrow::read_parquet(in_path)
log_msg(paste0("Loaded parquet: n=", nrow(df), " rows; p=", ncol(df), " columns."))

required <- c("motility_ratio","species","site","season","culture","experiment",
              "ntu","lux_exposure","cu_ug_l","toxin_exposure")
missing_req <- setdiff(required, names(df))
if (length(missing_req) > 0) stop("Missing required column(s): ", paste(missing_req, collapse = ", "))

if (!("video_file" %in% names(df))) {
  df <- df %>% mutate(video_file = paste0("row_", row_number()))
  log_msg("NOTE: video_file not present; created synthetic row id 'row_<n>'.")
}

# ---- Option A: experiment_id (same as Script 04) ----
df <- df %>%
  mutate(
    experiment_chr = as.character(experiment),
    culture_chr    = as.character(culture),
    experiment_id  = if_else(is.na(experiment_chr) | experiment_chr == "",
                             paste0("CTRL_", culture_chr),
                             experiment_chr)
  ) %>%
  mutate(
    species        = factor(species),
    site           = factor(site),
    season         = factor(season),
    culture        = factor(culture),
    experiment_id  = factor(experiment_id),
    toxin_exposure = factor(toxin_exposure)
  )

# ---- Beta squeeze ----
squeeze_beta <- function(y) {
  n <- sum(!is.na(y))
  if (n == 0) return(y)
  (y * (n - 1) + 0.5) / n
}
df <- df %>% mutate(motility_beta = squeeze_beta(motility_ratio))

# ---- Stabilising transforms (raw retained) ----
df <- df %>%
  mutate(
    cu_log10 = log10(cu_ug_l + 1),
    cu_z     = as.numeric(scale(cu_log10)),
    lux_z    = as.numeric(scale(lux_exposure)),
    ntu_z    = as.numeric(scale(ntu))
  )

log_msg("Created transformed/scaled covariates: cu_log10, cu_z, lux_z, ntu_z (raw columns retained).")

# ---- Design snapshot ----
design_ss <- df %>%
  count(site, species, name = "n_rows") %>%
  tidyr::complete(site, species, fill = list(n_rows = 0)) %>%
  arrange(site, species)

missing_cells <- design_ss %>% filter(n_rows == 0)
log_msg(paste0("Design site×species missing cells: ", nrow(missing_cells)))
if (nrow(missing_cells) > 0) {
  log_msg(paste0("Missing cells: ", paste(paste0(missing_cells$site, "×", missing_cells$species), collapse = "; ")))
}

# ---- Helpers ----
prepare_model_data <- function(formula, data, model_name) {
  vars <- unique(all.vars(formula))
  needed <- intersect(vars, names(data))
  
  df_sub <- data[, c("video_file", needed), drop = FALSE]
  cc <- stats::complete.cases(df_sub)
  
  dropped <- data %>%
    mutate(.complete = cc) %>%
    filter(!.complete) %>%
    select(video_file, all_of(needed))
  
  used <- data %>% filter(cc)
  
  out_drop <- here::here("outputs", paste0("model_rows_dropped_sens_", model_name, ".csv"))
  readr::write_csv(dropped, out_drop)
  
  log_msg(paste0(
    "Sensitivity model '", model_name, "': complete cases = ", nrow(used),
    " / ", nrow(data),
    " (dropped ", nrow(dropped), "). Wrote: ", out_drop
  ))
  
  list(data_used = used, dropped = dropped, dropped_path = out_drop)
}

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
      control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
    ),
    warning = w_handler
  )
  
  pdHess <- isTRUE(fit$sdr$pdHess)
  conv   <- fit$fit$convergence
  
  log_msg(paste0("  -> pdHess: ", pdHess, " ; convergence code: ", conv))
  if (length(warn) > 0) {
    log_msg(paste0("  -> warnings (", length(warn), "):"))
    for (i in seq_along(warn)) log_msg(paste0("     [", i, "] ", warn[i]))
  }
  
  list(fit = fit, pdHess = pdHess, conv = conv, warnings = warn)
}

fit_and_record <- function(name, formula, data) {
  log_msg(paste0("Fitting sensitivity model: ", name))
  log_msg(paste0("Formula: ", deparse(formula)))
  
  prep <- prepare_model_data(formula, data, model_name = name)
  res  <- fit_glmmTMB_logged(formula, prep$data_used, model_name = name)
  fit  <- res$fit
  
  # Save model
  rds_path <- here::here("models", paste0("ms4_motility_sens_", name, ".rds"))
  saveRDS(fit, rds_path)
  log_msg(paste0("  -> Saved model: ", rds_path))
  
  # Save summary
  sum_path <- here::here("outputs", paste0("model_summary_sens_", name, ".txt"))
  writeLines(capture.output(summary(fit)), sum_path)
  log_msg(paste0("  -> Wrote summary: ", sum_path))
  
  # DHARMa diagnostics
  diag_path <- here::here("outputs", paste0("model_diagnostics_sens_", name, ".txt"))
  sim <- DHARMa::simulateResiduals(fit, n = 500)
  diag_lines <- c(
    paste0("Sensitivity model diagnostics: ", name),
    paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Model:",
    deparse(formula),
    "",
    paste0("pdHess: ", res$pdHess),
    paste0("optimizer convergence code: ", res$conv),
    "",
    "DHARMa tests:"
  )
  diag_lines <- c(diag_lines, capture.output(print(DHARMa::testDispersion(sim))))
  diag_lines <- c(diag_lines, capture.output(print(DHARMa::testZeroInflation(sim))))
  diag_lines <- c(diag_lines, capture.output(print(DHARMa::testUniformity(sim))))
  writeLines(diag_lines, diag_path)
  log_msg(paste0("  -> Wrote diagnostics: ", diag_path))
  
  meta <- tibble::tibble(
    model_name = name,
    n_rows     = nrow(prep$data_used),
    n_culture  = dplyr::n_distinct(prep$data_used$culture),
    n_expid    = dplyr::n_distinct(prep$data_used$experiment_id),
    pdHess     = res$pdHess,
    conv_code  = res$conv,
    AIC        = AIC(fit),
    BIC        = BIC(fit),
    logLik     = as.numeric(logLik(fit))
  )
  
  list(fit = fit, meta = meta)
}

# ================================================================================
# Sensitivity model set
# ================================================================================
# SENS-1: Remove rank-deficient interaction (additive fixed effects)
sens1_formula <- motility_beta ~
  species + site + season +
  lux_z + cu_z + toxin_exposure +
  (1 | culture / experiment_id)

# SENS-2: Keep interaction but split random effects (often more stable)
sens2_formula <- motility_beta ~
  species * site + season +
  lux_z + cu_z + toxin_exposure +
  (1 | culture) + (1 | experiment_id)

# SENS-3: Additive fixed effects + split random effects (stable fallback)
sens3_formula <- motility_beta ~
  species + site + season +
  lux_z + cu_z + toxin_exposure +
  (1 | culture) + (1 | experiment_id)

# SENS-4: If toxin_exposure is the primary “SPM mechanism” driver, check a model
# where site/species are treated additively and toxin is left to do the work.
# (This is NOT for final inference, but helps diagnose whether interaction terms
# are causing instability.)
sens4_formula <- motility_beta ~
  species + site + season +
  lux_z + cu_z + toxin_exposure +
  (1 | culture) + (1 | experiment_id)

# (sens4 is same fixed structure as sens3 here — kept explicitly so you can later
# swap in alternative terms without altering the rest of the pipeline.)

# ---- Fit all sensitivity models ----
results <- list()
results[["additive_nestedRE"]]   <- fit_and_record("additive_nestedRE",   sens1_formula, df)
results[["interaction_splitRE"]] <- fit_and_record("interaction_splitRE", sens2_formula, df)
results[["additive_splitRE"]]    <- fit_and_record("additive_splitRE",    sens3_formula, df)
results[["additive_splitRE_2"]]  <- fit_and_record("additive_splitRE_2",  sens4_formula, df)

# ---- Comparison table ----
comparison <- dplyr::bind_rows(lapply(results, function(x) x$meta)) %>%
  arrange(AIC)

out_comp <- here::here("outputs", "sensitivity_model_comparison.csv")
readr::write_csv(comparison, out_comp)
log_msg(paste0("Wrote sensitivity model comparison: ", out_comp))

# ---- Fixed effects extracts (all sensitivity models) ----
fix_tables <- lapply(names(results), function(nm) {
  fit <- results[[nm]]$fit
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    mutate(model_name = nm)
}) %>% dplyr::bind_rows()

out_fix <- here::here("outputs", "sensitivity_fixed_effects_all.csv")
readr::write_csv(fix_tables, out_fix)
log_msg(paste0("Wrote combined sensitivity fixed effects: ", out_fix))

# ---- Human-readable summary ----
summary_path <- here::here("outputs", "sensitivity_summaries.txt")
lines <- c(
  "Manuscript 4 — Sensitivity models summary (05_models_sensitivity_motility.R)",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Goal: evaluate robustness to rank-deficiency and random-effects specification.",
  "",
  "Design notes:",
  paste0("- Missing site×species cells: ", nrow(missing_cells)),
  if (nrow(missing_cells) > 0) paste0("  * ", paste(paste0(missing_cells$site, "×", missing_cells$species), collapse = "; ")) else "  * None",
  "",
  "Model comparison (sorted by AIC):"
)
lines <- c(lines, capture.output(print(comparison, n = Inf, width = Inf)))
lines <- c(lines, "", "Interpretation guidance:",
           "- Prefer models with pdHess=TRUE and convergence code 0 for uncertainty reporting.",
           "- If key fixed-effect directions are stable across specs, inference is robust.",
           "- Primary inference remains Script 04; these are support checks.",
           "",
           "End of sensitivity summary.")
writeLines(lines, summary_path)
log_msg(paste0("Wrote sensitivity summary text: ", summary_path))

log_msg("05_models_sensitivity_motility.R complete.")


m <- readRDS("models/ms4_motility_core_glmmTMB.rds")
isTRUE(m$sdr$pdHess)
m$fit$convergence
