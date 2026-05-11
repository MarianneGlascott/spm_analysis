# =========================================================
# Script title: 17_format_supplementary_tables_word_01.R
# Project: SPM Analysis
# Purpose: Format Day 4 supplementary statistical tables
#          into manuscript-ready Word tables for A4 output.
#
# Inputs:
# - outputs/tables/day4/Supp_Table_S1_day4_candidate_model_comparison.csv
# - outputs/tables/day4/Supp_Table_S2_day4_preferred_model_diagnostics.csv
#
# Output:
# - outputs/tables/day4/Supplementary_day4_statistical_tables.docx
#
# Notes:
# - Tables are formatted in Times New Roman, 12 pt.
# - A4 landscape layout is used to improve readability.
# - Full audit-trail CSVs remain unchanged.
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 17: FORMAT SUPPLEMENTARY TABLES FOR WORD\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")

# ---------------------------------------------------------
# 1. Check setup
# ---------------------------------------------------------

required_objects <- c("project_root")

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Missing setup object(s):\n- ",
      paste(missing_objects, collapse = "\n- "),
      "\nPlease run 01_setup_packages_and_paths.R first."
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------
# 2. Package checks
# ---------------------------------------------------------

required_packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "officer",
  "flextable"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing package(s):\n- ",
      paste(missing_packages, collapse = "\n- "),
      "\nPlease install these before running this script."
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------
# 3. Paths
# ---------------------------------------------------------

dir_day4_tables <- file.path(project_root, "outputs", "tables", "day4")

file_s1 <- file.path(
  dir_day4_tables,
  "Supp_Table_S1_day4_candidate_model_comparison.csv"
)

file_s2 <- file.path(
  dir_day4_tables,
  "Supp_Table_S2_day4_preferred_model_diagnostics.csv"
)

file_s1_caption <- file.path(
  dir_day4_tables,
  "Supp_Table_S1_day4_candidate_model_comparison_caption.txt"
)

file_s2_caption <- file.path(
  dir_day4_tables,
  "Supp_Table_S2_day4_preferred_model_diagnostics_caption.txt"
)

file_docx <- file.path(
  dir_day4_tables,
  "Supplementary_day4_statistical_tables.docx"
)

# ---------------------------------------------------------
# 4. Helper functions
# ---------------------------------------------------------

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      paste0("Required file not found for ", label, ":\n", path),
      call. = FALSE
    )
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_caption <- function(path, fallback) {
  if (file.exists(path)) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  } else {
    fallback
  }
}

format_delta <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x == 0 ~ "0.00",
    TRUE ~ sprintf("%.2f", x)
  )
}

format_aic <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    TRUE ~ sprintf("%.2f", x)
  )
}

clean_model_name <- function(x) {
  x |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_replace_all("betabinomial", "beta-binomial") |>
    stringr::str_to_sentence()
}

make_manuscript_flextable <- function(df) {
  ft <- flextable::flextable(df)
  
  ft <- ft |>
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 12, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::align(align = "center", part = "header") |>
    flextable::align(align = "left", part = "body") |>
    flextable::valign(valign = "top", part = "all") |>
    flextable::padding(padding = 3, part = "all") |>
    flextable::border_remove() |>
    flextable::hline_top(
      border = officer::fp_border(color = "black", width = 1),
      part = "all"
    ) |>
    flextable::hline_bottom(
      border = officer::fp_border(color = "black", width = 1),
      part = "header"
    ) |>
    flextable::hline_bottom(
      border = officer::fp_border(color = "black", width = 1),
      part = "body"
    ) |>
    flextable::autofit()
  
  ft
}

# ---------------------------------------------------------
# 5. Read and prepare tables
# ---------------------------------------------------------

s1_raw <- read_required_csv(file_s1, "Supplementary Table S1")
s2_raw <- read_required_csv(file_s2, "Supplementary Table S2")

caption_s1 <- read_caption(
  file_s1_caption,
  "Supplementary Table S1. Candidate model comparison for Day 4 zoospore motility analyses."
)

caption_s2 <- read_caption(
  file_s2_caption,
  "Supplementary Table S2. Preferred Day 4 model summary and simulated residual diagnostics."
)

# ---------------------------------------------------------
# 6. Table S1: Candidate model comparison
# ---------------------------------------------------------

s1_table <- s1_raw |>
  dplyr::mutate(
    Experiment = dplyr::case_when(
      experiment == "Experiment 1" ~ "Exp. 1",
      experiment == "Experiment 2" ~ "Exp. 2",
      experiment == "Experiment 3" ~ "Exp. 3",
      experiment == "Experiment 4" ~ "Exp. 4",
      TRUE ~ experiment
    ),
    `Experimental pathway` = experiment_label,
    Model = clean_model_name(model),
    Family = dplyr::case_when(
      family == "betabinomial" ~ "Beta-binomial",
      family == "binomial" ~ "Binomial",
      TRUE ~ family
    ),
    df = as.character(df),
    AIC = format_aic(AIC),
    `Delta AIC` = format_delta(delta_aic),
    Selected = selected_model
  ) |>
  dplyr::select(
    Experiment,
    `Experimental pathway`,
    Model,
    Family,
    df,
    AIC,
    `Delta AIC`,
    Selected
  )

# ---------------------------------------------------------
# 7. Table S2: Preferred model diagnostics
# ---------------------------------------------------------

s2_table <- s2_raw |>
  dplyr::mutate(
    Experiment = dplyr::case_when(
      experiment == "Experiment 1" ~ "Exp. 1",
      experiment == "Experiment 2" ~ "Exp. 2",
      experiment == "Experiment 3" ~ "Exp. 3",
      experiment == "Experiment 4" ~ "Exp. 4",
      TRUE ~ experiment
    ),
    `Experimental pathway` = experiment_label,
    `Preferred model` = clean_model_name(preferred_model),
    Family = dplyr::case_when(
      preferred_family == "betabinomial" ~ "Beta-binomial",
      preferred_family == "binomial" ~ "Binomial",
      TRUE ~ preferred_family
    ),
    AIC = sprintf("%.2f", preferred_aic),
    `Delta AIC to next model` = sprintf("%.2f", next_best_delta_aic),
    `Model support` = model_selection_certainty,
    `Uniformity p` = uniformity_p,
    `Dispersion p` = dispersion_p,
    `Outlier p` = outlier_p,
    `Diagnostic interpretation` = diagnostic_interpretation
  ) |>
  dplyr::select(
    Experiment,
    `Experimental pathway`,
    `Preferred model`,
    Family,
    AIC,
    `Delta AIC to next model`,
    `Model support`,
    `Uniformity p`,
    `Dispersion p`,
    `Outlier p`,
    `Diagnostic interpretation`
  )

# ---------------------------------------------------------
# 8. Create flextables
# ---------------------------------------------------------

ft_s1 <- make_manuscript_flextable(s1_table)

# Narrower first columns / wider model column
ft_s1 <- ft_s1 |>
  flextable::width(j = "Experiment", width = 0.55) |>
  flextable::width(j = "Experimental pathway", width = 1.65) |>
  flextable::width(j = "Model", width = 1.85) |>
  flextable::width(j = "Family", width = 1.05) |>
  flextable::width(j = "df", width = 0.40) |>
  flextable::width(j = "AIC", width = 0.65) |>
  flextable::width(j = "Delta AIC", width = 0.75) |>
  flextable::width(j = "Selected", width = 0.70)

ft_s2 <- make_manuscript_flextable(s2_table)

ft_s2 <- ft_s2 |>
  flextable::width(j = "Experiment", width = 0.55) |>
  flextable::width(j = "Experimental pathway", width = 1.50) |>
  flextable::width(j = "Preferred model", width = 1.55) |>
  flextable::width(j = "Family", width = 1.05) |>
  flextable::width(j = "AIC", width = 0.60) |>
  flextable::width(j = "Delta AIC to next model", width = 0.90) |>
  flextable::width(j = "Model support", width = 1.00) |>
  flextable::width(j = "Uniformity p", width = 0.75) |>
  flextable::width(j = "Dispersion p", width = 0.75) |>
  flextable::width(j = "Outlier p", width = 0.70) |>
  flextable::width(j = "Diagnostic interpretation", width = 2.20)

# ---------------------------------------------------------
# 9. Build Word document
# ---------------------------------------------------------

doc <- officer::read_docx()

# A4 landscape section
landscape_section <- officer::prop_section(
  page_size = officer::page_size(orient = "landscape"),
  page_margins = officer::page_mar(
    top = 0.5,
    bottom = 0.5,
    left = 0.5,
    right = 0.5
  )
)

normal_text <- officer::fp_text(
  font.size = 12,
  font.family = "Times New Roman"
)

caption_text <- officer::fp_text(
  font.size = 12,
  font.family = "Times New Roman",
  bold = TRUE
)

doc <- doc |>
  officer::body_add_par(
    "Supplementary statistical tables",
    style = "heading 1"
  ) |>
  officer::body_add_par(
    "Formatted for A4 landscape layout, 12 pt Times New Roman.",
    style = "Normal"
  ) |>
  officer::body_end_section_continuous(value = landscape_section)

# Table S1
doc <- doc |>
  officer::body_add_fpar(
    officer::fpar(
      officer::ftext("Supplementary Table S1. ", prop = caption_text),
      officer::ftext(
        stringr::str_remove(caption_s1, "^Supplementary Table S1\\.\\s*"),
        prop = normal_text
      )
    )
  ) |>
  officer::body_add_flextable(ft_s1) |>
  officer::body_add_par("", style = "Normal") |>
  officer::body_add_break()

# Table S2
doc <- doc |>
  officer::body_add_fpar(
    officer::fpar(
      officer::ftext("Supplementary Table S2. ", prop = caption_text),
      officer::ftext(
        stringr::str_remove(caption_s2, "^Supplementary Table S2\\.\\s*"),
        prop = normal_text
      )
    )
  ) |>
  officer::body_add_flextable(ft_s2) |>
  officer::body_add_par("", style = "Normal")

print(doc, target = file_docx)

cat("\nWord document written to:\n")
cat(file_docx, "\n\n")

cat("SCRIPT 17 COMPLETE\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")