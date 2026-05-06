# =========================================================
# Script title: helpers_save_figures.R
# Project: SPM Analysis
# Author: Marianne Glascott
# Affiliation: School of Life Sciences, University of Sussex
# Manuscript: Manuscript 4
# Purpose: Provide central helper functions for saving
#          publication-ready figures in consistent formats,
#          dimensions, and filenames across the project.
# Inputs: ggplot objects or grid-compatible plots
# Outputs: Figure files written to outputs/figures/
#          in PDF, PNG, and TIFF formats
# Author: Marianne Glascott
# Date created: 24 February 2026
# Last updated: 24 March 2026
# Notes/dependencies:
# - Source via 01_setup_packages_and_paths.R
# - Designed to align with project output rules:
#   PDF, PNG, TIFF; 600 dpi for raster; stable filenames;
#   one- and two-column widths where relevant.
# - Uses dimensions defined in helpers_theme.R where
#   available; otherwise falls back to sensible defaults.
# =========================================================

message("Loading helper script: R/helpers_save_figures.R")

# ---------------------------------------------------------
# 1. Package checks
# ---------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required by helpers_save_figures.R.", call. = FALSE)
}

# ragg is optional but preferred for raster output quality
.use_ragg <- requireNamespace("ragg", quietly = TRUE)

# ---------------------------------------------------------
# 2. Dimension defaults
# ---------------------------------------------------------

# Use shared theme constants if already available; otherwise
# fall back to local defaults so this helper can still work.
if (!exists("fig_width_one_col", inherits = TRUE)) {
  fig_width_one_col <- 85 / 25.4
}

if (!exists("fig_width_two_col", inherits = TRUE)) {
  fig_width_two_col <- 178 / 25.4
}

if (!exists("fig_height_std", inherits = TRUE)) {
  fig_height_std <- 110 / 25.4
}

if (!exists("fig_height_tall", inherits = TRUE)) {
  fig_height_tall <- 140 / 25.4
}

if (!exists("fig_height_short", inherits = TRUE)) {
  fig_height_short <- 90 / 25.4
}

if (!exists("fig_dpi", inherits = TRUE)) {
  fig_dpi <- 600
}

# ---------------------------------------------------------
# 3. Output directory helper
# ---------------------------------------------------------

get_figure_dir <- function(subdir = NULL) {
  if (exists("dir_figures", inherits = TRUE)) {
    base_dir <- get("dir_figures", inherits = TRUE)
  } else if (requireNamespace("here", quietly = TRUE)) {
    base_dir <- here::here("outputs", "figures")
  } else {
    base_dir <- file.path("outputs", "figures")
  }

  if (!is.null(subdir) && nzchar(subdir)) {
    base_dir <- file.path(base_dir, subdir)
  }

  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }

  return(base_dir)
}

# ---------------------------------------------------------
# 4. Filename sanitising helper
# ---------------------------------------------------------

sanitize_filename <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---------------------------------------------------------
# 5. Width selector helper
# ---------------------------------------------------------

resolve_figure_width <- function(width = c("one_col", "two_col", "custom"),
                                 custom_width = NULL) {
  width <- match.arg(width)

  if (width == "one_col") {
    return(fig_width_one_col)
  }

  if (width == "two_col") {
    return(fig_width_two_col)
  }

  if (width == "custom") {
    if (is.null(custom_width) || !is.numeric(custom_width) || length(custom_width) != 1) {
      stop("For width = 'custom', provide a single numeric custom_width in inches.", call. = FALSE)
    }
    return(custom_width)
  }

  stop("Invalid width specification.", call. = FALSE)
}

resolve_figure_height <- function(height = c("standard", "tall", "short", "custom"),
                                  custom_height = NULL) {
  height <- match.arg(height)

  if (height == "standard") {
    return(fig_height_std)
  }

  if (height == "tall") {
    return(fig_height_tall)
  }

  if (height == "short") {
    return(fig_height_short)
  }

  if (height == "custom") {
    if (is.null(custom_height) || !is.numeric(custom_height) || length(custom_height) != 1) {
      stop("For height = 'custom', provide a single numeric custom_height in inches.", call. = FALSE)
    }
    return(custom_height)
  }

  stop("Invalid height specification.", call. = FALSE)
}

# ---------------------------------------------------------
# 6. Core device writers
# ---------------------------------------------------------

save_plot_pdf <- function(plot, filename, width, height, bg = "white") {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = bg
  )
}

save_plot_png <- function(plot, filename, width, height, dpi = fig_dpi, bg = "white") {
  if (.use_ragg) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_png,
      bg = bg
    )
  } else {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "png",
      bg = bg
    )
  }
}

save_plot_tiff <- function(plot, filename, width, height, dpi = fig_dpi, bg = "white") {
  if (.use_ragg) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw",
      bg = bg
    )
  } else {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "tiff",
      compression = "lzw",
      bg = bg
    )
  }
}

# ---------------------------------------------------------
# 7. Main figure-saving helper
# ---------------------------------------------------------

save_figure_set <- function(plot,
                            figure_name,
                            subdir = NULL,
                            width = c("one_col", "two_col", "custom"),
                            height = c("standard", "tall", "short", "custom"),
                            custom_width = NULL,
                            custom_height = NULL,
                            dpi = fig_dpi,
                            save_pdf = TRUE,
                            save_png = TRUE,
                            save_tiff = TRUE,
                            bg = "white",
                            quiet = FALSE) {
  if (missing(plot) || is.null(plot)) {
    stop("plot must be supplied.", call. = FALSE)
  }

  if (missing(figure_name) || !nzchar(figure_name)) {
    stop("figure_name must be a non-empty string.", call. = FALSE)
  }

  width_value <- resolve_figure_width(
    width = match.arg(width),
    custom_width = custom_width
  )

  height_value <- resolve_figure_height(
    height = match.arg(height),
    custom_height = custom_height
  )

  out_dir <- get_figure_dir(subdir = subdir)
  safe_name <- sanitize_filename(figure_name)

  file_pdf  <- file.path(out_dir, paste0(safe_name, ".pdf"))
  file_png  <- file.path(out_dir, paste0(safe_name, ".png"))
  file_tiff <- file.path(out_dir, paste0(safe_name, ".tiff"))

  if (save_pdf) {
    save_plot_pdf(
      plot = plot,
      filename = file_pdf,
      width = width_value,
      height = height_value,
      bg = bg
    )
  }

  if (save_png) {
    save_plot_png(
      plot = plot,
      filename = file_png,
      width = width_value,
      height = height_value,
      dpi = dpi,
      bg = bg
    )
  }

  if (save_tiff) {
    save_plot_tiff(
      plot = plot,
      filename = file_tiff,
      width = width_value,
      height = height_value,
      dpi = dpi,
      bg = bg
    )
  }

  saved_files <- c(
    if (save_pdf) file_pdf,
    if (save_png) file_png,
    if (save_tiff) file_tiff
  )

  if (!quiet) {
    message("Saved figure set for: ", safe_name)
    for (f in saved_files) {
      message("  - ", f)
    }
  }

  invisible(saved_files)
}

# ---------------------------------------------------------
# 8. Dual-width helper
# ---------------------------------------------------------

save_figure_both_widths <- function(plot,
                                    figure_name,
                                    subdir = NULL,
                                    height = c("standard", "tall", "short", "custom"),
                                    custom_height = NULL,
                                    dpi = fig_dpi,
                                    save_pdf = TRUE,
                                    save_png = TRUE,
                                    save_tiff = TRUE,
                                    bg = "white",
                                    quiet = FALSE) {
  if (missing(figure_name) || !nzchar(figure_name)) {
    stop("figure_name must be a non-empty string.", call. = FALSE)
  }

  files_one <- save_figure_set(
    plot = plot,
    figure_name = paste0(figure_name, "_one_col"),
    subdir = subdir,
    width = "one_col",
    height = match.arg(height),
    custom_height = custom_height,
    dpi = dpi,
    save_pdf = save_pdf,
    save_png = save_png,
    save_tiff = save_tiff,
    bg = bg,
    quiet = quiet
  )

  files_two <- save_figure_set(
    plot = plot,
    figure_name = paste0(figure_name, "_two_col"),
    subdir = subdir,
    width = "two_col",
    height = match.arg(height),
    custom_height = custom_height,
    dpi = dpi,
    save_pdf = save_pdf,
    save_png = save_png,
    save_tiff = save_tiff,
    bg = bg,
    quiet = quiet
  )

  invisible(c(files_one, files_two))
}

# ---------------------------------------------------------
# 9. Figure caption helper
# ---------------------------------------------------------

write_figure_caption_md <- function(figure_name,
                                    caption_text,
                                    subdir = NULL,
                                    overwrite = TRUE) {
  if (missing(figure_name) || !nzchar(figure_name)) {
    stop("figure_name must be provided.", call. = FALSE)
  }

  if (missing(caption_text) || !nzchar(caption_text)) {
    stop("caption_text must be provided.", call. = FALSE)
  }

  out_dir <- get_figure_dir(subdir = subdir)
  safe_name <- sanitize_filename(figure_name)
  file_md <- file.path(out_dir, paste0(safe_name, "_caption.md"))

  if (file.exists(file_md) && !overwrite) {
    stop("Caption file already exists and overwrite = FALSE.", call. = FALSE)
  }

  writeLines(caption_text, con = file_md)

  message("Saved figure caption markdown: ", file_md)

  invisible(file_md)
}

# ---------------------------------------------------------
# 10. Save figure + caption wrapper
# ---------------------------------------------------------

save_figure_with_caption <- function(plot,
                                     figure_name,
                                     caption_text = NULL,
                                     subdir = NULL,
                                     width = c("one_col", "two_col", "custom"),
                                     height = c("standard", "tall", "short", "custom"),
                                     custom_width = NULL,
                                     custom_height = NULL,
                                     dpi = fig_dpi,
                                     save_pdf = TRUE,
                                     save_png = TRUE,
                                     save_tiff = TRUE,
                                     bg = "white",
                                     quiet = FALSE,
                                     write_caption = TRUE) {
  saved_files <- save_figure_set(
    plot = plot,
    figure_name = figure_name,
    subdir = subdir,
    width = match.arg(width),
    height = match.arg(height),
    custom_width = custom_width,
    custom_height = custom_height,
    dpi = dpi,
    save_pdf = save_pdf,
    save_png = save_png,
    save_tiff = save_tiff,
    bg = bg,
    quiet = quiet
  )

  caption_file <- NULL

  if (write_caption && !is.null(caption_text) && nzchar(caption_text)) {
    caption_file <- write_figure_caption_md(
      figure_name = figure_name,
      caption_text = caption_text,
      subdir = subdir
    )
  }

  invisible(list(
    figure_files = saved_files,
    caption_file = caption_file
  ))
}

# ---------------------------------------------------------
# 11. Simple log helper
# ---------------------------------------------------------

log_saved_figure <- function(figure_name,
                             files_saved,
                             log_path = NULL) {
  if (is.null(log_path)) {
    if (exists("dir_logs", inherits = TRUE)) {
      log_path <- file.path(
        get("dir_logs", inherits = TRUE),
        "figure_save_log.txt"
      )
    } else if (requireNamespace("here", quietly = TRUE)) {
      log_path <- here::here("outputs", "logs", "figure_save_log.txt")
    } else {
      log_path <- file.path("outputs", "logs", "figure_save_log.txt")
    }
  }

  log_dir <- dirname(log_path)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  lines <- c(
    paste0("[", timestamp, "] ", figure_name),
    paste0("  ", files_saved),
    ""
  )

  cat(lines, file = log_path, sep = "\n", append = TRUE)

  invisible(log_path)
}

# ---------------------------------------------------------
# 12. Convenience wrapper with logging
# ---------------------------------------------------------

save_and_log_figure <- function(plot,
                                figure_name,
                                caption_text = NULL,
                                subdir = NULL,
                                width = c("one_col", "two_col", "custom"),
                                height = c("standard", "tall", "short", "custom"),
                                custom_width = NULL,
                                custom_height = NULL,
                                dpi = fig_dpi,
                                save_pdf = TRUE,
                                save_png = TRUE,
                                save_tiff = TRUE,
                                bg = "white",
                                quiet = FALSE,
                                write_caption = TRUE,
                                log_path = NULL) {
  result <- save_figure_with_caption(
    plot = plot,
    figure_name = figure_name,
    caption_text = caption_text,
    subdir = subdir,
    width = match.arg(width),
    height = match.arg(height),
    custom_width = custom_width,
    custom_height = custom_height,
    dpi = dpi,
    save_pdf = save_pdf,
    save_png = save_png,
    save_tiff = save_tiff,
    bg = bg,
    quiet = quiet,
    write_caption = write_caption
  )

  log_saved_figure(
    figure_name = figure_name,
    files_saved = result$figure_files,
    log_path = log_path
  )

  invisible(result)
}

message("helpers_save_figures.R loaded successfully.")
message("Available helpers: save_figure_set(), save_figure_both_widths(), save_figure_with_caption(), save_and_log_figure()")