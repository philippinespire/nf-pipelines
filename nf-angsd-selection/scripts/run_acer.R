#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT NAME:  run_acer.R
# PURPOSE:      Executes an iterative genome-wide selection scan and temporal Ne 
#               estimation using chi-squared and CMH tests from the ACER package 
#               across historical and modern population timepoints.
#
# USAGE:        Rscript run_acer.R --hist_mafs <hist1.mafs.gz,hist2.mafs.gz> \
#                                  --mod_mafs <mod1.mafs.gz,mod2.mafs.gz> \
#                                  --region_names <reg1,reg2> [options]
#
# REQUIRED ARGUMENTS:
#   --hist_mafs     Comma-separated list of historical ANGSD .mafs.gz file paths
#   --mod_mafs      Comma-separated list of modern ANGSD .mafs.gz file paths
#   --region_names  Comma-separated list of region/population names corresponding
#                   to the order of the MAF files
#
# OPTIONAL ARGUMENTS:
#   --out_dir       Output results directory (Default: ./results/selection)
#   --helpers       Path to acer_helpers.R script (Default: scripts/acer_helpers.R)
#   --generations   Elapsed generations between sampling points (Default: 114)
#   --fdr_cutoff    False discovery rate threshold for candidate loci (Default: 0.05)
#   --max_rounds    Maximum iterations for Ne/selection scan convergence (Default: 20)
#   --n_boot        Bootstrap iterations for Ne confidence intervals (Default: 1000)
#   --min_ind       Minimum individual count threshold per site (Default: 4)
#
# OUTPUTS:       tsv files for iteration summary, final test results, CMH 
#               statistics, and Ne bootstrap confidence intervals.
# ==============================================================================

suppressPackageStartupMessages({
  library(ACER)
  library(boot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  hist_mafs      = NULL,
  mod_mafs       = NULL,
  region_names   = NULL,
  out_dir        = "./results/selection",
  helpers        = "scripts/acer_helpers.R",
  generations    = 114,
  fdr_cutoff     = 0.05,
  max_rounds     = 20,
  n_boot         = 1000,
  min_ind        = 4
)

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

if (is.null(opt$hist_mafs) || is.null(opt$mod_mafs) || is.null(opt$region_names)) {
  stop("Error: Must provide --hist_mafs, --mod_mafs, and --region_names.", call.=FALSE)
}

if (!file.exists(opt$helpers)) stop(paste("Helper script not found at:", opt$helpers))
source(opt$helpers)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Automatically derive ne_generations and test_generations
ne_generations   <- opt$generations
test_generations <- c(0, opt$generations - 1)

hist_files <- strsplit(opt$hist_mafs, ",")[[1]]
mod_files  <- strsplit(opt$mod_mafs, ",")[[1]]
reg_names  <- strsplit(opt$region_names, ",")[[1]]

merged_regions_list <- list()
regions <- list()

for (i in seq_along(reg_names)) {
  r_name <- reg_names[i]
  hist_df <- read.table(gzfile(hist_files[i]), header = TRUE, stringsAsFactors = FALSE)
  mod_df  <- read.table(gzfile(mod_files[i]), header = TRUE, stringsAsFactors = FALSE)
  
  af_col_hist <- ifelse("knownEM" %in% names(hist_df), "knownEM", "freq")
  af_col_mod  <- ifelse("knownEM" %in% names(mod_df), "knownEM", "freq")
  
  merged <- merge(hist_df, mod_df, by = c("chromo", "position"), suffixes = c("_H", "_M"))
  
  # Ensure site allele compatibility
  valid_sites <- (merged$major_H == merged$major_M & merged$minor_H == merged$minor_M) | 
                 (merged$major_H == merged$minor_M & merged$minor_H == merged$major_M)
  merged <- merged[valid_sites, , drop = FALSE]
  
  is_flipped <- merged$major_H == merged$minor_M & merged$minor_H == merged$major_M
  merged$AF_mod_adj <- merged[[paste0(af_col_mod, "_M")]]
  merged$AF_mod_adj[is_flipped] <- 1 - merged$AF_mod_adj[is_flipped]
  
  reg_df <- data.frame(CHR = merged$chromo, BP = merged$position)
  reg_df[[paste0(r_name, "_A_AF")]] <- merged[[paste0(af_col_hist, "_H")]]
  reg_df[[paste0(r_name, "_C_AF")]] <- merged$AF_mod_adj
  reg_df[[paste0(r_name, "_A_N")]]  <- merged$nInd_H
  reg_df[[paste0(r_name, "_C_N")]]  <- merged$nInd_M
  
  merged_regions_list[[i]] <- reg_df
  regions[[r_name]] <- c(A = paste0(r_name, "_A"), C = paste0(r_name, "_C"))
}

# Merge and filter across regions
df_unfiltered <- merged_regions_list[[1]]
if (length(merged_regions_list) > 1) {
  for (i in 2:length(merged_regions_list)) {
    df_unfiltered <- merge(df_unfiltered, merged_regions_list[[i]], by = c("CHR", "BP"), all = FALSE)
  }
}

af_cols <- grep("_AF$", names(df_unfiltered), value = TRUE)
df_unfiltered <- df_unfiltered[complete.cases(df_unfiltered[, af_cols]), , drop = FALSE]

all_n_cols <- grep("_N$", names(df_unfiltered), value = TRUE)
df <- df_unfiltered[apply(df_unfiltered[, all_n_cols, drop = FALSE] >= opt$min_ind, 1, all), , drop = FALSE]

# Ensure polymorphic sites only (drop sites fixed at 0 or 1 across all samples)
poly_mask <- apply(df[, af_cols, drop = FALSE], 1, function(x) !all(x == 0) && !all(x == 1))
df <- df[poly_mask, , drop = FALSE]

message(sprintf("Total polymorphic SNPs passing filters: %d", nrow(df)))

# Iterative selection scan
iterative_result <- run_iterative_selection(
  df = df, regions = regions, ne_generations = ne_generations,
  test_generations = test_generations, fdr_cutoff = opt$fdr_cutoff, max_rounds = opt$max_rounds
)

iteration_summary <- build_iteration_summary(iterative_result, regions)
final_outputs <- build_final_outputs(iterative_result$final, regions)

write.table(iteration_summary, file.path(opt$out_dir, "iteration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_outputs$test_results, file.path(opt$out_dir, "final_test_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_outputs$cmh_results_full, file.path(opt$out_dir, "final_cmh_results_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Bootstrap Ne
final_selected_idx <- iterative_result$final$selected_idx
df_neutral_final <- if (length(final_selected_idx) > 0) df[-final_selected_idx, , drop = FALSE] else df
neutral_data_final <- build_region_data(df_neutral_final, regions)
ne_boot <- bootstrap_ne_by_region(neutral_data_final$mafs, neutral_data_final$covs, regions, generations = ne_generations, n_boot = opt$n_boot)
write.table(ne_boot, file.path(opt$out_dir, "ne_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)