# Project:     SPM Analysis - Impact on kelp zoospore motility
# Script:      08_NTU as continuous (Log scale).R
# Author:      Marianne Glascott
# Date:        2026-02-13
# Experiment   8.2
# Purpose:
#   Test hypothesis: H1 turbidity
# Input:
#   data_raw/All_data_cleaned.csv
# Outputs (written to disk; required):
#
# ==============================================================================
#===============================================================================
# 0. Print helper function
#===============================================================================
save_pub_figure <- function(plot, filename_base,
                            width = 6, height = 5,
                            dpi = 600) {
  
  base_path <- here::here("figures_pub", filename_base)
  
  # PDF (vector)
  ggsave(
    paste0(base_path, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  # PNG (high resolution)
  ggsave(
    paste0(base_path, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  # TIFF (journal submission standard)
  ggsave(
    paste0(base_path, ".tif"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}


#===============================================================================
library(dplyr)
library(ggplot2)
library(arrow)

df_p4 <- arrow::read_parquet("data_processed/ms4_p4_analysis.parquet")

df_8 <- df_p4 %>%
  filter(experiment_chr == "8.2")

range(df_8$ntu, na.rm = TRUE)
#===============================================================================
# 1. Load + sanity + subset
#===============================================================================
library(dplyr)
library(ggplot2)
library(readr)
library(arrow)
library(here)
library(scales)

# Optional: your project styling helpers (if you have them)
# source(here::here("R", "kelp_style.R"))
# source(here::here("R", "helpers.R"))

df_p4 <- arrow::read_parquet(here::here("data_processed", "ms4_p4_analysis.parquet"))

df_8 <- df_p4 %>%
  filter(experiment_chr == "8.2") %>%
  mutate(
    toxin_exposure = as.character(toxin_exposure),
    ntu = as.numeric(ntu),
    motility_ratio = as.numeric(motility_ratio),
    mobile_cell_count = as.numeric(mobile_cell_count),
    total_cells = as.numeric(total_cells)
  )

stopifnot(all(c("CONTROL","SPM") %in% unique(df_8$toxin_exposure)))

df_8 %>% count(toxin_exposure)
range(df_8$ntu, na.rm = TRUE)
summary(df_8$motility_ratio)
#===============================================================================
# 2. Wilson CI helper + summary table
#===============================================================================
wilson_ci <- function(x, n, conf = 0.95) {
  # Wilson score interval for a proportion
  z <- qnorm(1 - (1 - conf)/2)
  p <- x / n
  denom <- 1 + (z^2)/n
  center <- (p + (z^2)/(2*n)) / denom
  half <- (z * sqrt((p*(1-p) + (z^2)/(4*n))/n)) / denom
  tibble(lwr = pmax(0, center - half), upr = pmin(1, center + half))
}

# Summarise by NTU level (pooling replicates using counts)
sum_8 <- df_8 %>%
  group_by(toxin_exposure, ntu) %>%
  summarise(
    n_reps = n(),
    mobile = sum(mobile_cell_count, na.rm = TRUE),
    total  = sum(total_cells, na.rm = TRUE),
    motility = ifelse(total > 0, mobile / total, NA_real_),
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$mobile, .$total)) %>%
  rename(ci_low = lwr, ci_high = upr)

sum_8
#===============================================================================
# 3. Plot (raw points + summary + CI ribbon)
#===============================================================================
# A simple, consistent x-transform (handles ntu = 0)
df_8 <- df_8 %>% mutate(ntu_plot = log10(ntu + 1))
sum_8 <- sum_8 %>% mutate(ntu_plot = log10(ntu + 1))

# If you have a species palette already, we can drop it in.
# For now we keep it clean and use default ggplot colours or a manual palette later.

p1 <- ggplot() +
  # Raw replicate points (transparent, jittered slightly)
  geom_point(
    data = df_8,
    aes(x = ntu_plot, y = motility_ratio, shape = toxin_exposure),
    alpha = 0.5,
    position = position_jitter(width = 0.02, height = 0)
  ) +
  # Summary CI ribbon
  geom_ribbon(
    data = sum_8,
    aes(x = ntu_plot, ymin = ci_low, ymax = ci_high, group = toxin_exposure),
    alpha = 0.20
  ) +
  # Summary line + points
  geom_line(
    data = sum_8,
    aes(x = ntu_plot, y = motility, group = toxin_exposure),
    linewidth = 0.8
  ) +
  geom_point(
    data = sum_8,
    aes(x = ntu_plot, y = motility),
    size = 2
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(
    name = "Turbidity (NTU; log10(NTU + 1))",
    breaks = log10(c(0, 100, 500, 1000, 2000, 4000) + 1),
    labels = c("0", "100", "500", "1000", "2000", "4000")
  ) +
  labs(
    y = "Motile fraction",
    shape = "Treatment"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold")
  )

p1
#===============================================================================
# 4. Save outputs (PNG + PDF)
#===============================================================================
dir.create(here::here("figures_pub"), showWarnings = FALSE)

ggsave(here::here("figures_pub", "Fig_NTU_gradient_exp8_2.pdf"), p1,
       width = 170, height = 115, units = "mm", dpi = 300)

ggsave(here::here("figures_pub", "Fig_NTU_gradient_exp8_2.png"), p1,
       width = 170, height = 115, units = "mm", dpi = 300)
#===============================================================================
# 5. Descriptive statistics
#===============================================================================
# ---- Descriptive statistics: Exp 8.2 NTU gradient ----

out_path_summary <- here::here("outputs", "exp8_ntugradient_motility_summary.txt")
dir.create(dirname(out_path_summary), showWarnings = FALSE)

summary_lines <- c(
  "Experiment 8.2 — NTU Gradient Motility Summary",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Rows: ", nrow(df_8)),
  "",
  "Overall motility_ratio summary:",
  capture.output(summary(df_8$motility_ratio)),
  "",
  "By toxin_exposure:",
  capture.output(
    df_8 %>%
      group_by(toxin_exposure) %>%
      summarise(
        n = n(),
        min = min(motility_ratio, na.rm = TRUE),
        q1  = quantile(motility_ratio, 0.25, na.rm = TRUE),
        median = median(motility_ratio, na.rm = TRUE),
        mean = mean(motility_ratio, na.rm = TRUE),
        q3  = quantile(motility_ratio, 0.75, na.rm = TRUE),
        max = max(motility_ratio, na.rm = TRUE)
      )
  )
)

writeLines(summary_lines, out_path_summary)

message("Saved summary to: ", out_path_summary)

# manuscript friendly
readr::write_csv(
  df_8 %>%
    group_by(toxin_exposure, ntu) %>%
    summarise(
      n = n(),
      mean_motility = mean(motility_ratio, na.rm = TRUE),
      sd_motility   = sd(motility_ratio, na.rm = TRUE),
      .groups = "drop"
    ),
  here::here("outputs", "exp8_ntu_descriptives_by_level.csv")
)

#===============================================================================
# 6. Continuous or categorical?
#===============================================================================
cor.test(df_8$ntu, df_8$motility_ratio, method = "spearman")

# Spearman correlation
spearman_test <- cor.test(
  df_8$ntu,
  df_8$motility_ratio,
  method = "spearman"
)
out_path <- here::here("outputs", "exp8_spearman_ntu_vs_motility.txt")

writeLines(
  c(
    "Experiment 8.2 — Spearman correlation",
    paste0("Timestamp: ", Sys.time()),
    "",
    capture.output(print(spearman_test))
  ),
  con = out_path
)

# saved output as .csv
spearman_summary <- tibble::tibble(
  method = spearman_test$method,
  rho = unname(spearman_test$estimate),
  statistic_S = unname(spearman_test$statistic),
  p_value = spearman_test$p.value,
  conf_low = spearman_test$conf.int[1],
  conf_high = spearman_test$conf.int[2],
  n = sum(complete.cases(df_8$ntu, df_8$motility_ratio))
)

readr::write_csv(
  spearman_summary,
  here::here("outputs", "exp8_spearman_summary.csv")
)


#===============================================================================
# 7. Fit linear model
#===============================================================================
lm_8 <- lm(motility_ratio ~ ntu, data = df_8)

summary_lm_8 <- summary(lm_8)

out_lm_path <- here::here("outputs", "exp8_lm_motility_vs_ntu.txt")

writeLines(
  c(
    "Experiment 8.2 — Linear model: motility_ratio ~ ntu",
    paste0("Timestamp: ", Sys.time()),
    "",
    capture.output(print(summary_lm_8))
  ),
  out_lm_path
)

lm_summary_df <- tibble::tibble(
  term = rownames(coef(summary_lm_8)),
  estimate = coef(summary_lm_8)[, "Estimate"],
  std_error = coef(summary_lm_8)[, "Std. Error"],
  t_value = coef(summary_lm_8)[, "t value"],
  p_value = coef(summary_lm_8)[, "Pr(>|t|)"]
)

readr::write_csv(
  lm_summary_df,
  here::here("outputs", "exp8_lm_summary.csv")
)

#===============================================================================
# 8. Diagnostic plots
#===============================================================================
png(here::here("outputs", "exp8_lm_diagnostics.png"), width = 1200, height = 1000, res = 150)
par(mfrow = c(2, 2))
plot(lm_8)
dev.off()

#===============================================================================
# 9. Create clean plot
#===============================================================================
library(ggplot2)

df_8 <- df_8 %>%
  mutate(log10_ntu = log10(ntu + 1))

p_exp8 <- ggplot(df_8, aes(x = log10_ntu, y = motility_ratio)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", colour = "black", fill = "grey70") +
  labs(
    x = "Turbidity (log10 NTU + 1)",
    y = "Motile fraction"
  ) +
  theme_classic(base_size = 14)

save_pub_figure(
  plot = p_exp8,
  filename_base = "Fig1_exp8_motility_vs_turbidity",
  width = 6,
  height = 5,
  dpi = 600
)
 # confirm save

file.exists("figures_pub/Fig1_exp8_motility_vs_turbidity.tif")
file.info("figures_pub/Fig1_exp8_motility_vs_turbidity.tif")$size
file.exists("figures_pub/Fig1_exp8_motility_vs_turbidity.pdf")
file.info("figures_pub/Fig1_exp8_motility_vs_turbidity.pdf")$size
file.exists("figures_pub/Fig1_exp8_motility_vs_turbidity.png")
file.info("figures_pub/Fig1_exp8_motility_vs_turbidity.png")$size

#===============================================================================
# 10. Figure export
#===============================================================================
save_pub_figure <- function(plot, filename_base,
                            width = 6, height = 5,
                            dpi = 600) {
  
  base_path <- here::here("figures_pub", filename_base)
  
  # PDF (vector)
  ggsave(
    paste0(base_path, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  # PNG (high resolution)
  ggsave(
    paste0(base_path, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  # TIFF (journal submission standard)
  ggsave(
    paste0(base_path, ".tif"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

save_pub_figure(
  plot = p_exp8,
  filename_base = "Fig1_exp8_motility_vs_turbidity",
  width = 6,
  height = 5,
  dpi = 600
)

#===============================================================================
# 11. Figure refinement
#===============================================================================
kelp_colour <- "#c09c0e"

p_exp8 <- ggplot(df_8, aes(x = ntu, y = motility_ratio)) +
  
  # Raw data
  geom_point(
    colour = kelp_colour,
    alpha = 0.6,
    size = 2
  ) +
  
  # Smoothed trend (loess for visual only)
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = kelp_colour,
    fill = "grey80",
    linewidth = 1
  ) +
  
  scale_x_continuous(
    "Turbidity (NTU)",
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  scale_y_continuous(
    "Motile fraction",
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA)
  )
# add spearman result
p_exp8 <- p_exp8 +
  annotate(
    "text",
    x = max(df_8$ntu) * 0.6,
    y = 0.95,
    label = "Spearman ρ = -0.006\np = 0.96",
    hjust = 0,
    size = 3.5
  )

#===============================================================================
# 11. Figure refinement
#===============================================================================
kelp_colour <- "#c09c0e"

p_exp8 <- ggplot(df_8, aes(x = ntu, y = motility_ratio)) +
  
  # Raw data
  geom_point(
    colour = kelp_colour,
    alpha = 0.6,
    size = 2
  ) +
  
  # Smoothed trend (loess for visual only)
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = kelp_colour,
    fill = "grey80",
    linewidth = 1
  ) +
  
  scale_x_continuous(
    "Turbidity (NTU)",
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  scale_y_continuous(
    "Motile fraction",
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA)
  )
# add spearman result
p_exp8 <- p_exp8 +
  annotate(
    "text",
    x = max(df_8$ntu) * 0.6,
    y = 0.95,
    label = "Spearman ρ = -0.006\np = 0.96",
    hjust = 0,
    size = 3.5
  )

print(p_exp8)

#===============================================================================
# 12. Save figure - Exp 8 NTU vs Motility
#===============================================================================
#===============================================================================
# 12. Save Figure – Exp 8 NTU vs Motility
#===============================================================================

out_dir <- here::here("figures_pub")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

base_name <- "Fig1_Exp8_NTU_vs_Motility"

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".png")),
  plot = p_exp8,
  width = 160,
  height = 120,
  units = "mm",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".pdf")),
  plot = p_exp8,
  width = 160,
  height = 120,
  units = "mm"
)

ggsave(
  filename = file.path(out_dir, paste0(base_name, ".tif")),
  plot = p_exp8,
  width = 160,
  height = 120,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

cat("Saved figure to:", out_dir, "\n")
