# =========================================================
# SCRIPT 12: SUPPLEMENTARY OUTPUTS (REVISED)
# Focus: supporting robustness + transparency
# =========================================================

cat("\n========================================================\n")
cat("SCRIPT 12: SUPPLEMENTARY OUTPUTS\n")
cat("========================================================\n\n")

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# --- Load model comparisons ---
mod_exp4 <- read_if_exists("outputs/tables/09_exp4_field_spm_model_comparison.csv")
mod_exp4_day4 <- read_if_exists("outputs/models/exp4_field_spm/exp4_field_spm_day4_model_comparison.csv")

mod_exp4$dataset <- "Full"
mod_exp4_day4$dataset <- "Day 4"

combined <- dplyr::bind_rows(mod_exp4, mod_exp4_day4)

# --- Plot ---
p <- ggplot2::ggplot(combined) +
  ggplot2::geom_col(
    ggplot2::aes(x = model, y = delta_aic, fill = dataset),
    position = "dodge"
  ) +
  ggplot2::labs(
    title = "Experiment 4 model comparison (sensitivity analysis)",
    y = "ΔAIC",
    x = "Model"
  )

ggplot2::ggsave("outputs/figures/supplementary/FigS3_model_comparison.png", p, width = 7, height = 5)

cat("Supplementary Fig S3 saved.\n")