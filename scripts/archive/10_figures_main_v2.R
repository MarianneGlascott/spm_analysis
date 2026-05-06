# =========================================================
# SCRIPT 10: FIGURES MAIN (REVISED FOR MANUSCRIPT 4)
# Focus: clear hypothesis-driven figures, with emphasis on H4
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 10: FIGURES MAIN\n")
cat("========================================================\n\n")

# --- Load data ---
dat <- readRDS(file.path(dir_data_derived, "ms4_analysis_derived.rds"))

# Restrict to MS4 experiments
dat <- dat |>
  dplyr::filter(experiment_num %in% c(8.2, 9.2, 10.2, 11.2))

# --- Load predictions ---
read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

pred_exp1 <- read_if_exists(file.path(project_root,"outputs/models/exp1_light/exp1_light_preferred_model_predictions.csv"))
pred_exp2 <- read_if_exists(file.path(project_root,"outputs/models/exp2_defined_particles/exp2_defined_particles_preferred_model_predictions.csv"))
pred_exp3 <- read_if_exists(file.path(project_root,"outputs/models/exp3_brake_size/exp3_brake_size_preferred_model_predictions.csv"))
pred_exp4 <- read_if_exists(file.path(project_root,"outputs/models/exp4_field_spm/exp4_field_spm_preferred_model_predictions.csv"))
pred_exp4_day4 <- read_if_exists(file.path(project_root,"outputs/models/exp4_field_spm/exp4_field_spm_day4_preferred_model_predictions.csv"))

# --- Common prep ---
prep_dat <- function(df) {
  df |>
    dplyr::mutate(
      total_cells = mobile_cell_count + stationary_cell_count,
      motility_ratio = mobile_cell_count / total_cells
    )
}

# --- EXP 4 (PRIMARY FIGURE) ---
exp4 <- dat |> dplyr::filter(experiment_num == 8.2) |> prep_dat()

p_full <- ggplot2::ggplot(exp4) +
  ggplot2::geom_point(
    ggplot2::aes(log10(ntu + 1), motility_ratio),
    alpha = 0.25
  ) +
  ggplot2::geom_line(
    data = pred_exp4,
    ggplot2::aes(log10_ntu_plus1, fit_response),
    colour = "#0072B2", linewidth = 1
  ) +
  ggplot2::facet_wrap(~days_from_start) +
  ggplot2::labs(
    subtitle = "Full dataset (Days 4–12)",
    x = "log10(NTU + 1)",
    y = "Motile fraction"
  )

p_day4 <- ggplot2::ggplot(exp4 |> dplyr::filter(days_from_start == 4)) +
  ggplot2::geom_point(
    ggplot2::aes(log10(ntu + 1), motility_ratio),
    alpha = 0.25
  ) +
  ggplot2::geom_line(
    data = pred_exp4_day4,
    ggplot2::aes(log10_ntu_plus1, fit_response),
    colour = "#0072B2", linewidth = 1
  ) +
  ggplot2::labs(
    subtitle = "Day 4 sensitivity",
    x = "log10(NTU + 1)",
    y = "Motile fraction"
  )

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_fig5 <- p_full / p_day4 +
    patchwork::plot_annotation(
      title = "Field-derived SPM gradient",
      subtitle = "Nonlinear NTU-response and robustness to Day 4 restriction"
    )
} else {
  p_fig5 <- p_full
}

ggplot2::ggsave("outputs/figures/main/Fig5_field_spm_gradient.png", p_fig5, width = 8, height = 10, dpi = 600)

cat("Fig 5 saved.\n")