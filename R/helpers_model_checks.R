# =========================================================
# Script title: helpers_model_checks.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Define reusable helper functions for fitting,
#          checking, summarising, comparing, predicting from,
#          and saving models across the analysis pipeline.
# Inputs: Analysis-ready data frames and fitted model objects
# Outputs: Reusable helper functions; model objects, tables,
#          logs, and diagnostics written by downstream scripts
# Date created: 24 February 2026
# Last updated: 24 March 2026
# Notes/dependencies:
# - Source via 01_setup_packages_and_paths.R
# - Intended primarily for glmmTMB binomial models fitted to
#   cbind(mobile_cell_count, stationary_cell_count)
# - Designed to support model summaries, diagnostics,
#   predictions, model comparisons, and reproducible saving
# =========================================================

message("Loading helper script: R/helpers_model_checks.R")

# ---------------------------------------------------------
# 1. Package checks
# ---------------------------------------------------------

required_model_packages <- c(
  "glmmTMB",
  "broom.mixed",
  "dplyr",
  "readr",
  "tibble"
)

missing_model_packages <- required_model_packages[
  !vapply(required_model_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_model_packages) > 0) {
  stop(
    paste0(
      "helpers_model_checks.R requires the following package(s): ",
      paste(missing_model_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

.has_dharma <- requireNamespace("DHARMa", quietly = TRUE)
.has_performance <- requireNamespace("performance", quietly = TRUE)
.has_emmeans <- requireNamespace("emmeans", quietly = TRUE)

# ---------------------------------------------------------
# 2. Output directory helper
# ---------------------------------------------------------

get_model_dir <- function(subdir = NULL) {
  if (exists("dir_models", inherits = TRUE)) {
    base_dir <- get("dir_models", inherits = TRUE)
  } else if (requireNamespace("here", quietly = TRUE)) {
    base_dir <- here::here("outputs", "models")
  } else {
    base_dir <- file.path("outputs", "models")
  }

  if (!is.null(subdir) && nzchar(subdir)) {
    base_dir <- file.path(base_dir, subdir)
  }

  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }

  return(base_dir)
}

get_model_log_dir <- function() {
  if (exists("dir_logs", inherits = TRUE)) {
    base_dir <- get("dir_logs", inherits = TRUE)
  } else if (requireNamespace("here", quietly = TRUE)) {
    base_dir <- here::here("outputs", "logs")
  } else {
    base_dir <- file.path("outputs", "logs")
  }

  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }

  return(base_dir)
}

# ---------------------------------------------------------
# 3. Filename sanitising helper
# ---------------------------------------------------------

sanitize_model_filename <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---------------------------------------------------------
# 4. Small messaging helper
# ---------------------------------------------------------

model_message <- function(...) {
  message("[model] ", paste0(..., collapse = ""))
}

# ---------------------------------------------------------
# 5. Formula / response checks
# ---------------------------------------------------------

assert_primary_count_response <- function(data,
                                          mobile_col = "mobile_cell_count",
                                          stationary_col = "stationary_cell_count") {
  required_cols <- c(mobile_col, stationary_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing primary response column(s): ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_no_missing_model_vars <- function(data, vars) {
  vars <- unique(vars)
  vars <- vars[vars %in% names(data)]

  missing_counts <- vapply(vars, function(v) sum(is.na(data[[v]])), integer(1))
  bad <- missing_counts[missing_counts > 0]

  if (length(bad) > 0) {
    warning(
      paste0(
        "Missing values remain in model variable(s): ",
        paste(names(bad), bad, sep = "=", collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------
# 6. Core model fitting wrappers
# ---------------------------------------------------------

fit_glmmtmb_binomial <- function(formula,
                                 data,
                                 ziformula = ~0,
                                 dispformula = ~1,
                                 weights = NULL,
                                 REML = FALSE,
                                 save_name = NULL,
                                 save_subdir = NULL,
                                 quiet = FALSE) {
  if (missing(formula)) {
    stop("A model formula must be provided.", call. = FALSE)
  }

  if (missing(data) || !is.data.frame(data)) {
    stop("A data frame must be provided.", call. = FALSE)
  }

  assert_primary_count_response(data)

  if (!quiet) {
    model_message("Fitting glmmTMB binomial model...")
  }

  fit <- glmmTMB::glmmTMB(
    formula = formula,
    data = data,
    family = stats::binomial(link = "logit"),
    ziformula = ziformula,
    dispformula = dispformula,
    weights = weights,
    REML = REML
  )

  if (!quiet) {
    model_message("Model fit complete.")
  }

  if (!is.null(save_name)) {
    save_model_object(
      model = fit,
      model_name = save_name,
      subdir = save_subdir,
      quiet = quiet
    )
  }

  fit
}

# ---------------------------------------------------------
# 7. Save model object helper
# ---------------------------------------------------------

save_model_object <- function(model,
                              model_name,
                              subdir = NULL,
                              quiet = FALSE) {
  if (missing(model) || is.null(model)) {
    stop("A fitted model must be supplied.", call. = FALSE)
  }

  if (missing(model_name) || !nzchar(model_name)) {
    stop("model_name must be provided.", call. = FALSE)
  }

  out_dir <- get_model_dir(subdir = subdir)
  safe_name <- sanitize_model_filename(model_name)
  out_file <- file.path(out_dir, paste0(safe_name, ".rds"))

  saveRDS(model, out_file)

  if (!quiet) {
    model_message("Saved model object: ", out_file)
  }

  invisible(out_file)
}

# ---------------------------------------------------------
# 8. Tidy fixed-effect summary helper
# ---------------------------------------------------------

tidy_model_fixed <- function(model,
                             conf.int = TRUE,
                             conf.level = 0.95,
                             exponentiate = FALSE) {
  out <- broom.mixed::tidy(
    model,
    effects = "fixed",
    conf.int = conf.int,
    conf.level = conf.level,
    exponentiate = exponentiate
  )

  out <- out |>
    dplyr::rename(
      std_error = std.error,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value
    ) |>
    dplyr::mutate(
      effect_scale = ifelse(exponentiate, "odds_ratio", "logit"),
      term = as.character(term)
    )

  out
}

# ---------------------------------------------------------
# 9. Random effects summary helper
# ---------------------------------------------------------

tidy_model_random <- function(model) {
  out <- broom.mixed::tidy(
    model,
    effects = "ran_pars",
    conf.int = FALSE
  )

  out <- out |>
    dplyr::rename(
      std_error = std.error
    )

  out
}

# ---------------------------------------------------------
# 10. Model metadata helper
# ---------------------------------------------------------

extract_model_metadata <- function(model,
                                   model_name = NULL) {
  tibble::tibble(
    model_name = ifelse(is.null(model_name), NA_character_, model_name),
    class = paste(class(model), collapse = "; "),
    family = tryCatch(model$family$family, error = function(e) NA_character_),
    link = tryCatch(model$family$link, error = function(e) NA_character_),
    nobs = tryCatch(stats::nobs(model), error = function(e) NA_integer_),
    logLik = tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_),
    AIC = tryCatch(stats::AIC(model), error = function(e) NA_real_),
    BIC = tryCatch(stats::BIC(model), error = function(e) NA_real_)
  )
}

# ---------------------------------------------------------
# 11. Fixed + metadata summary table
# ---------------------------------------------------------

build_model_term_summary <- function(model,
                                     model_name = NULL,
                                     conf.level = 0.95,
                                     exponentiate = FALSE) {
  fixed_terms <- tidy_model_fixed(
    model = model,
    conf.int = TRUE,
    conf.level = conf.level,
    exponentiate = exponentiate
  )

  meta <- extract_model_metadata(
    model = model,
    model_name = model_name
  )

  fixed_terms |>
    dplyr::mutate(
      model = meta$model_name[[1]],
      family = meta$family[[1]],
      link = meta$link[[1]]
    )
}

# ---------------------------------------------------------
# 12. Diagnostics helpers
# ---------------------------------------------------------

run_dharma_diagnostics <- function(model,
                                   n_sim = 1000,
                                   quiet = FALSE) {
  if (!.has_dharma) {
    warning("Package 'DHARMa' not available; DHARMa diagnostics skipped.", call. = FALSE)
    return(NULL)
  }

  if (!quiet) {
    model_message("Running DHARMa simulation diagnostics...")
  }

  sim_res <- DHARMa::simulateResiduals(
    fittedModel = model,
    n = n_sim
  )

  out <- list(
    simulation = sim_res,
    uniformity = tryCatch(DHARMa::testUniformity(sim_res), error = function(e) e),
    dispersion = tryCatch(DHARMa::testDispersion(sim_res), error = function(e) e),
    outliers = tryCatch(DHARMa::testOutliers(sim_res), error = function(e) e),
    zero_inflation = tryCatch(DHARMa::testZeroInflation(sim_res), error = function(e) e)
  )

  if (!quiet) {
    model_message("DHARMa diagnostics complete.")
  }

  out
}

run_performance_checks <- function(model,
                                   quiet = FALSE) {
  if (!.has_performance) {
    warning("Package 'performance' not available; performance checks skipped.", call. = FALSE)
    return(NULL)
  }

  if (!quiet) {
    model_message("Running performance checks...")
  }

  out <- list(
    check_model = tryCatch(performance::check_model(model), error = function(e) e),
    check_overdispersion = tryCatch(performance::check_overdispersion(model), error = function(e) e),
    check_zeroinflation = tryCatch(performance::check_zeroinflation(model), error = function(e) e),
    r2 = tryCatch(performance::r2(model), error = function(e) e)
  )

  if (!quiet) {
    model_message("Performance checks complete.")
  }

  out
}

# ---------------------------------------------------------
# 13. Compact diagnostics summary table
# ---------------------------------------------------------

summarise_diagnostics <- function(dharma_checks = NULL,
                                  performance_checks = NULL) {
  rows <- list()

  if (!is.null(dharma_checks)) {
    rows <- c(
      rows,
      list(tibble::tibble(
        diagnostic = "DHARMa uniformity",
        result = paste(class(dharma_checks$uniformity)[1]),
        note = if (inherits(dharma_checks$uniformity, "htest")) dharma_checks$uniformity$method else "See object"
      )),
      list(tibble::tibble(
        diagnostic = "DHARMa dispersion",
        result = paste(class(dharma_checks$dispersion)[1]),
        note = if (inherits(dharma_checks$dispersion, "htest")) dharma_checks$dispersion$method else "See object"
      )),
      list(tibble::tibble(
        diagnostic = "DHARMa outliers",
        result = paste(class(dharma_checks$outliers)[1]),
        note = if (inherits(dharma_checks$outliers, "htest")) dharma_checks$outliers$method else "See object"
      )),
      list(tibble::tibble(
        diagnostic = "DHARMa zero inflation",
        result = paste(class(dharma_checks$zero_inflation)[1]),
        note = if (inherits(dharma_checks$zero_inflation, "htest")) dharma_checks$zero_inflation$method else "See object"
      ))
    )
  }

  if (!is.null(performance_checks)) {
    rows <- c(
      rows,
      list(tibble::tibble(
        diagnostic = "Performance overdispersion",
        result = paste(class(performance_checks$check_overdispersion)[1]),
        note = "See object"
      )),
      list(tibble::tibble(
        diagnostic = "Performance zero inflation",
        result = paste(class(performance_checks$check_zeroinflation)[1]),
        note = "See object"
      )),
      list(tibble::tibble(
        diagnostic = "Performance R2",
        result = paste(class(performance_checks$r2)[1]),
        note = "See object"
      ))
    )
  }

  if (length(rows) == 0) {
    return(tibble::tibble(
      diagnostic = character(),
      result = character(),
      note = character()
    ))
  }

  dplyr::bind_rows(rows)
}

# ---------------------------------------------------------
# 14. Model comparison helper
# ---------------------------------------------------------

compare_models_aic <- function(...) {
  models <- list(...)

  if (length(models) < 2) {
    stop("Provide at least two fitted models for comparison.", call. = FALSE)
  }

  model_names <- names(models)
  if (is.null(model_names) || any(!nzchar(model_names))) {
    model_names <- paste0("model_", seq_along(models))
  }

  out <- tibble::tibble(
    model = model_names,
    AIC = vapply(models, stats::AIC, numeric(1)),
    BIC = vapply(models, stats::BIC, numeric(1)),
    logLik = vapply(models, function(m) as.numeric(stats::logLik(m)), numeric(1)),
    nobs = vapply(models, stats::nobs, numeric(1))
  ) |>
    dplyr::arrange(AIC) |>
    dplyr::mutate(
      delta_aic = AIC - min(AIC, na.rm = TRUE),
      delta_bic = BIC - min(BIC, na.rm = TRUE)
    )

  out
}

likelihood_ratio_compare <- function(model_small,
                                     model_large) {
  stats::anova(model_small, model_large)
}

# ---------------------------------------------------------
# 15. Prediction helpers
# ---------------------------------------------------------

make_prediction_grid <- function(data,
                                 focal_terms,
                                 at = list(),
                                 n_numeric = 100) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (missing(focal_terms) || length(focal_terms) == 0) {
    stop("At least one focal term must be provided.", call. = FALSE)
  }

  focal_terms <- unique(focal_terms)

  grid_list <- lapply(focal_terms, function(v) {
    if (!v %in% names(data)) {
      stop(paste0("Focal term not found in data: ", v), call. = FALSE)
    }

    if (v %in% names(at)) {
      return(at[[v]])
    }

    x <- data[[v]]

    if (is.numeric(x)) {
      return(seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n_numeric))
    }

    unique(stats::na.omit(x))
  })

  names(grid_list) <- focal_terms

  newdata <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  remaining_vars <- setdiff(names(data), focal_terms)

  for (v in remaining_vars) {
    x <- data[[v]]

    if (v %in% names(at)) {
      newdata[[v]] <- at[[v]][1]
    } else if (is.factor(x)) {
      newdata[[v]] <- levels(x)[1]
    } else if (is.character(x)) {
      newdata[[v]] <- unique(stats::na.omit(x))[1]
    } else if (is.numeric(x)) {
      newdata[[v]] <- stats::median(x, na.rm = TRUE)
    } else {
      newdata[[v]] <- unique(stats::na.omit(x))[1]
    }
  }

  for (v in names(data)) {
    if (is.factor(data[[v]]) && v %in% names(newdata)) {
      newdata[[v]] <- factor(newdata[[v]], levels = levels(data[[v]]))
    }
  }

  newdata
}

predict_glmmtmb_response <- function(model,
                                     newdata,
                                     conf.level = 0.95,
                                     re.form = NA) {
  pred_link <- stats::predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    re.form = re.form
  )

  z <- stats::qnorm((1 + conf.level) / 2)

  out <- tibble::as_tibble(newdata) |>
    dplyr::mutate(
      fit_link = as.numeric(pred_link$fit),
      se_link = as.numeric(pred_link$se.fit),
      conf_low_link = fit_link - z * se_link,
      conf_high_link = fit_link + z * se_link,
      fit_response = stats::plogis(fit_link),
      conf_low_response = stats::plogis(conf_low_link),
      conf_high_response = stats::plogis(conf_high_link)
    )

  out
}

# ---------------------------------------------------------
# 16. emmeans helper
# ---------------------------------------------------------

get_emmeans_table <- function(model,
                              specs,
                              type = "response") {
  if (!.has_emmeans) {
    warning("Package 'emmeans' not available; emmeans skipped.", call. = FALSE)
    return(NULL)
  }

  emm <- emmeans::emmeans(model, specs = specs, type = type)
  out <- as.data.frame(emm)

  tibble::as_tibble(out)
}

# ---------------------------------------------------------
# 17. Save model summary tables
# ---------------------------------------------------------

write_model_summary_csv <- function(data,
                                    file_stem,
                                    subdir = NULL,
                                    quiet = FALSE) {
  out_dir <- get_model_dir(subdir = subdir)
  safe_name <- sanitize_model_filename(file_stem)
  out_file <- file.path(out_dir, paste0(safe_name, ".csv"))

  readr::write_csv(data, out_file)

  if (!quiet) {
    model_message("Saved model summary CSV: ", out_file)
  }

  invisible(out_file)
}

write_model_summary_rds <- function(object,
                                    file_stem,
                                    subdir = NULL,
                                    quiet = FALSE) {
  out_dir <- get_model_dir(subdir = subdir)
  safe_name <- sanitize_model_filename(file_stem)
  out_file <- file.path(out_dir, paste0(safe_name, ".rds"))

  saveRDS(object, out_file)

  if (!quiet) {
    model_message("Saved model summary RDS: ", out_file)
  }

  invisible(out_file)
}

# ---------------------------------------------------------
# 18. Save diagnostics bundle
# ---------------------------------------------------------

save_diagnostics_bundle <- function(diagnostics,
                                    bundle_name,
                                    subdir = NULL,
                                    quiet = FALSE) {
  out_dir <- get_model_dir(subdir = subdir)
  safe_name <- sanitize_model_filename(bundle_name)
  out_file <- file.path(out_dir, paste0(safe_name, ".rds"))

  saveRDS(diagnostics, out_file)

  if (!quiet) {
    model_message("Saved diagnostics bundle: ", out_file)
  }

  invisible(out_file)
}

# ---------------------------------------------------------
# 19. Log helper
# ---------------------------------------------------------

log_model_run <- function(model_name,
                          note_lines,
                          log_file = NULL) {
  if (is.null(log_file)) {
    log_file <- file.path(
      get_model_log_dir(),
      "model_run_log.txt"
    )
  }

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  lines <- c(
    paste0("[", timestamp, "] ", model_name),
    paste0("  ", note_lines),
    ""
  )

  cat(lines, file = log_file, sep = "\n", append = TRUE)

  invisible(log_file)
}

# ---------------------------------------------------------
# 20. One-stop model workflow helper
# ---------------------------------------------------------

fit_check_save_model <- function(formula,
                                 data,
                                 model_name,
                                 model_subdir = NULL,
                                 save_terms_csv = TRUE,
                                 save_meta_csv = TRUE,
                                 save_diagnostics = TRUE,
                                 run_dharma = TRUE,
                                 run_performance = TRUE,
                                 conf.level = 0.95,
                                 exponentiate = FALSE,
                                 quiet = FALSE) {
  fit <- fit_glmmtmb_binomial(
    formula = formula,
    data = data,
    save_name = model_name,
    save_subdir = model_subdir,
    quiet = quiet
  )

  term_summary <- build_model_term_summary(
    model = fit,
    model_name = model_name,
    conf.level = conf.level,
    exponentiate = exponentiate
  )

  meta_summary <- extract_model_metadata(
    model = fit,
    model_name = model_name
  )

  dharma_checks <- NULL
  performance_checks <- NULL

  if (run_dharma) {
    dharma_checks <- run_dharma_diagnostics(fit, quiet = quiet)
  }

  if (run_performance) {
    performance_checks <- run_performance_checks(fit, quiet = quiet)
  }

  diagnostics_summary <- summarise_diagnostics(
    dharma_checks = dharma_checks,
    performance_checks = performance_checks
  )

  if (save_terms_csv) {
    write_model_summary_csv(
      data = term_summary,
      file_stem = paste0(model_name, "_fixed_effects"),
      subdir = model_subdir,
      quiet = quiet
    )
  }

  if (save_meta_csv) {
    write_model_summary_csv(
      data = meta_summary,
      file_stem = paste0(model_name, "_metadata"),
      subdir = model_subdir,
      quiet = quiet
    )

    write_model_summary_csv(
      data = diagnostics_summary,
      file_stem = paste0(model_name, "_diagnostics_summary"),
      subdir = model_subdir,
      quiet = quiet
    )
  }

  if (save_diagnostics) {
    save_diagnostics_bundle(
      diagnostics = list(
        dharma = dharma_checks,
        performance = performance_checks,
        diagnostics_summary = diagnostics_summary
      ),
      bundle_name = paste0(model_name, "_diagnostics_bundle"),
      subdir = model_subdir,
      quiet = quiet
    )
  }

  log_model_run(
    model_name = model_name,
    note_lines = c(
      paste0("Formula: ", deparse(formula)),
      paste0("Rows used: ", tryCatch(stats::nobs(fit), error = function(e) NA_integer_)),
      paste0("AIC: ", round(tryCatch(stats::AIC(fit), error = function(e) NA_real_), 2)),
      paste0("BIC: ", round(tryCatch(stats::BIC(fit), error = function(e) NA_real_), 2))
    )
  )

  invisible(list(
    model = fit,
    term_summary = term_summary,
    metadata = meta_summary,
    dharma = dharma_checks,
    performance = performance_checks,
    diagnostics_summary = diagnostics_summary
  ))
}

# ---------------------------------------------------------
# 21. Convenience sensitivity helper
# ---------------------------------------------------------

fit_and_compare_models <- function(model_list,
                                   quiet = FALSE) {
  if (!is.list(model_list) || length(model_list) < 2) {
    stop("model_list must be a named list with at least two fitted models.", call. = FALSE)
  }

  if (is.null(names(model_list)) || any(!nzchar(names(model_list)))) {
    stop("model_list must be a named list of fitted models.", call. = FALSE)
  }

  comp <- compare_models_aic(!!!model_list)

  if (!quiet) {
    model_message("Model comparison complete.")
  }

  comp
}

# ---------------------------------------------------------
# 22. Load message
# ---------------------------------------------------------

message("helpers_model_checks.R loaded successfully.")
message("Available fit helpers: fit_glmmtmb_binomial(), fit_check_save_model()")
message("Available summary helpers: tidy_model_fixed(), tidy_model_random(), build_model_term_summary(), extract_model_metadata()")
message("Available diagnostic helpers: run_dharma_diagnostics(), run_performance_checks(), summarise_diagnostics()")
message("Available comparison/prediction helpers: compare_models_aic(), make_prediction_grid(), predict_glmmtmb_response(), get_emmeans_table()")

