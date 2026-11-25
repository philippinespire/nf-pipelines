#!/usr/bin/env Rscript

# DATE:         November 2026
# AUTHOR:       Malin Pinsky with assistance from Google Gemini
# ==============================================================================
# Script: plot_tp.R
# Description: Calculates mean nucleotide diversity (pi / tP per site) and 
#              Watterson's theta (theta_W / tW per site) with 95% bootstrap 
#              confidence intervals across historical and modern timepoints for 
#              a given region. Outputs summary TSV tables and comparison plots 
#              for both estimators.
# 
# Usage: 
#   Rscript plot_tp.R --historic_file <hist.pestPG> --modern_file <mod.pestPG> [options]
#
# Required Arguments:
#   --historic_file <path> : Path to the historic ANGSD .pestPG theta file
#   --modern_file <path>   : Path to the modern ANGSD .pestPG theta file
#
# Optional Arguments:
#   --region <string>        : Region or population label (Default: "Region")
#   --out_pi_png <path>      : Output path for pi PNG plot (Default: "pi_historic_vs_modern.png")
#   --out_pi_tsv <path>      : Output path for pi summary table (Default: "pi_summary.tsv")
#   --out_theta_png <path>   : Output path for Watterson's theta PNG plot 
#                              (Default: derived from --out_pi_png, e.g., "tW_historic_vs_modern.png")
#   --out_theta_tsv <path>   : Output path for Watterson's theta summary table 
#                              (Default: derived from --out_pi_tsv, e.g., "tW_summary.tsv")
#   --n_boot <int>           : Number of bootstrap iterations for CIs (Default: 10000)
#
# Example:
#   Rscript plot_tp.R \
#     --historic_file=Bali_historic.pestPG \
#     --modern_file=Bali_modern.pestPG \
#     --region=Bali \
#     --out_pi_png=pi_historic_vs_modern_Bali.png \
#     --out_pi_tsv=pi_summary_Bali.tsv \
#     --out_theta_png=tW_historic_vs_modern_Bali.png \
#     --out_theta_tsv=tW_summary_Bali.tsv
#
# Dependencies: Base R only (no external packages required)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  historic_file = NULL,
  modern_file   = NULL,
  region        = "Region",
  out_pi_png    = "pi_historic_vs_modern.png",
  out_pi_tsv    = "pi_summary.tsv",
  out_theta_png = "tW_historic_vs_modern.png",
  out_theta_tsv = "tW_summary.tsv",
  n_boot        = 10000
)

# Parse command line options
if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      parts <- if (grepl("=", arg)) strsplit(sub("^--", "", arg), "=")[[1]] else c(sub("^--", "", arg), args[i + 1])
      key <- parts[1]
      val <- parts[2]
      i <- if (grepl("=", arg)) i + 1 else i + 2
      if (key %in% names(opt)) {
        opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(val) else val
      }
    } else { i <- i + 1 }
  }
}

if (is.null(opt$historic_file) || is.null(opt$modern_file)) {
  stop("Error: Must provide both --historic_file and --modern_file.", call. = FALSE)
}

# Read in pi and Watterson's theta data from ANGSD .pestPG files
read_tp <- function(path, era_label) {
  dat <- read.table(
    file = path, header = TRUE, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE, comment.char = "", quote = ""
  )

  missing_cols <- setdiff(c("tP", "tW", "nSites"), names(dat))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s) [%s] in file: %s",
                 paste(missing_cols, collapse = ", "), path))
  }

  tp_vals     <- suppressWarnings(as.numeric(dat$tP))
  tw_vals     <- suppressWarnings(as.numeric(dat$tW))
  nsites_vals <- suppressWarnings(as.numeric(dat$nSites))
  keep        <- is.finite(tp_vals) & is.finite(tw_vals) & is.finite(nsites_vals) & nsites_vals > 0
  
  data.frame(
    era    = era_label,
    tP     = tp_vals[keep],
    tW     = tw_vals[keep],
    nSites = nsites_vals[keep],
    stringsAsFactors = FALSE
  )
}

bootstrap_theta_stat <- function(vals, nsites_vals, stat_col = "mean_tP", n_boot = 10000, seed = 42) {
  set.seed(seed)
  n <- length(vals)
  obs_val <- sum(vals) / sum(nsites_vals)
  
  boot_stats <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    sum(vals[idx]) / sum(nsites_vals[idx])
  }, numeric(1))
  
  out <- data.frame(
    n        = n,
    est      = obs_val,
    ci_lower = as.numeric(quantile(boot_stats, 0.025)),
    ci_upper = as.numeric(quantile(boot_stats, 0.975))
  )
  names(out)[2] <- stat_col
  out
}

summarize_metric <- function(df, value_col, mean_col, era_order, n_boot) {
  summary_list <- lapply(era_order, function(era_name) {
    rows <- df$era == era_name
    out  <- bootstrap_theta_stat(df[[value_col]][rows], df$nSites[rows], stat_col = mean_col, n_boot = n_boot)
    out$era <- era_name
    out
  })
  summary_df     <- do.call(rbind, summary_list)
  summary_df$era <- factor(summary_df$era, levels = era_order)
  summary_df[, c("era", "n", mean_col, "ci_lower", "ci_upper")]
}

plot_metric <- function(summary_df, mean_col, out_png, ylab_text, main_title) {
  era_order <- levels(summary_df$era)
  cols      <- c("Historical" = "#1f78b4", "Modern" = "#e31a1c")
  x_pos     <- seq_along(era_order)
  y_lim     <- range(summary_df$ci_lower, summary_df$ci_upper)
  y_pad     <- if (diff(y_lim) == 0) 0.01 else 0.08 * diff(y_lim)

  png(filename = out_png, width = 1800, height = 1200, res = 200, type = "cairo")
  par(mar = c(6, 6, 4, 2) + 0.1, las = 1)

  plot(
    x = x_pos,
    y = summary_df[[mean_col]],
    pch = 19, cex = 2,
    col = cols[as.character(summary_df$era)],
    xaxt = "n", xlab = "Era", ylab = ylab_text,
    ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)
  )

  arrows(
    x0 = x_pos, y0 = summary_df$ci_lower,
    x1 = x_pos, y1 = summary_df$ci_upper,
    angle = 90, code = 3, length = 0.08, lwd = 2,
    col = cols[as.character(summary_df$era)]
  )

  axis(1, at = x_pos, labels = era_order)
  grid(nx = NA, ny = NULL, col = "gray85", lty = "dotted")
  title(main = main_title)

  invisible(dev.off())
  message(sprintf("Plot saved to %s", out_png))
}

historic_df <- read_tp(opt$historic_file, "Historical")
modern_df   <- read_tp(opt$modern_file, "Modern")
all_df      <- rbind(historic_df, modern_df)
era_order <- c("Historical", "Modern")

# 1. Nucleotide Diversity (pi = tP / nSites)
pi_summary <- summarize_metric(all_df, "tP", "mean_tP", era_order, opt$n_boot)
write.table(pi_summary, file = opt$out_pi_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
plot_metric(
  pi_summary, "mean_tP", opt$out_pi_png,
  ylab_text  = "Mean pi (tP per site)",
  main_title = paste("Nucleotide Diversity (pi) -", opt$region)
)

# 2. Watterson's Theta (theta_W = tW / nSites)
theta_summary <- summarize_metric(all_df, "tW", "mean_tW", era_order, opt$n_boot)
write.table(theta_summary, file = opt$out_theta_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
plot_metric(
  theta_summary, "mean_tW", opt$out_theta_png,
  ylab_text  = "Mean Watterson's theta (tW per site)",
  main_title = paste("Watterson's Theta (theta_W) -", opt$region)
)