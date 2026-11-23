#!/usr/bin/env Rscript

# goal is to calculate Reference Bias (Ref / (Ref + Alt)), and Strand Bias (p < 0.05 fraction) across eras.
# also plot individual heterozygosity (Ho) vs. mean depth across callable sites for each individual.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# Force headless PNG rendering without X11
options(bitmapType = "cairo")

# read arguments
args <- commandArgs(trailingOnly = TRUE)
alist <- list()
i <- 1
while (i <= length(args)) {
  if (startsWith(args[i], "--")) {
    if (grepl("=", args[i])) {
      kv <- strsplit(sub("^--", "", args[i]), "=")[[1]]
      alist[[kv[1]]] <- paste(kv[-1], collapse = "=")
      i <- i + 1
    } else {
      key <- sub("^--", "", args[i])
      alist[[key]] <- args[i + 1]
      i <- i + 2
    }
  } else {
    i <- i + 1
  }
}
maf_files   <- unlist(strsplit(alist$mafs, ","))
count_files <- unlist(strsplit(alist$counts, ","))
if (is.null(alist$snpstat) || !nzchar(alist$snpstat)) {
  stop("Missing required --snpstat argument.")
}
snpstat_files <- unlist(strsplit(alist$snpstat, ","))

ind_sfs_files <- unlist(strsplit(alist$ind_sfs, ","))
meta <- fread(alist$samplesheet)
meta[, `:=`(
  sample = trimws(as.character(sample)),
  region = trimws(as.character(region)),
  era_norm = ifelse(
    grepl("^(historic|historical|hist|era1)$", trimws(tolower(as.character(era)))),
    "Historic", "Modern"
  )
)]
ind_depth_list <- list()

results_list <- list()
maf_dt_list  <- list()

# Read the target SNP list
snps <- fread(alist$snps, select = 1:2, col.names = c("chromo", "position"))

for (i in seq_along(maf_files)) {
  m_file <- maf_files[i]
  
  # Extract population and era from filename (e.g., Bali_historic.mafs.gz)
  base_name <- basename(m_file)
  pop_era   <- gsub("\\.mafs\\.gz$", "", base_name)
  target_era       <- ifelse(grepl("_(historic|historical|hist|era1)$", pop_era, ignore.case = TRUE), "Historic", "Modern")
  target_pop       <- sub("_(historic|historical|hist|era1|modern|mod|era2)$", "", pop_era, ignore.case = TRUE)
  
  # Match corresponding counts file by exact basename (Nextflow .collect() order is not guaranteed)
  c_match <- count_files[basename(count_files) == paste0(pop_era, ".counts.gz")]
  if (length(c_match) != 1 || !file.exists(c_match)) {
    stop(sprintf("Could not uniquely match counts file for %s", pop_era))
  }
  c_file <- c_match[1]

  # Per-individual mean depth across all callable sites in <pop>_<era>.counts.gz
  depth_cmd <- sprintf(
    "zcat %s | awk 'NR>1 {n_sites++; for(j=1;j<=NF;j++) col[j]+=$j} END {n_ind=NF/4; for(k=0;k<n_ind;k++) {d=col[4*k+1]+col[4*k+2]+col[4*k+3]+col[4*k+4]; printf \"%%d\\t%%.6f\\n\", k, (n_sites>0 ? d/n_sites : 0)}}'",
    shQuote(c_file)
  )
  cohort_depths <- fread(cmd = depth_cmd, col.names = c("ind_idx", "mean_depth"))
  cohort_meta <- meta[tolower(region) == tolower(target_pop) & era_norm == target_era]

  if (nrow(cohort_meta) != nrow(cohort_depths)) {
  stop(sprintf("Sample count mismatch for %s_%s: samplesheet has %d, counts.gz has %d",
              target_pop, target_era, nrow(cohort_meta), nrow(cohort_depths)))
  }
  cohort_depths[, `:=`(sample = cohort_meta$sample, pop = target_pop, era = target_era)]
  ind_depth_list[[i]] <- cohort_depths[, .(sample, pop, era, mean_depth)]

  # Read MAF data
  maf <- fread(cmd = paste("zcat", shQuote(m_file)), select = c("chromo", "position", "major", "minor", "ref"))

  # Read per-individual base counts via awk, skipping line 1 header (ind0_A ind0_C ind0_G ind0_T ...)
  awk_cmd <- sprintf(
    "zcat %s | awk 'NR>1 {a=0;c=0;g=0;t=0; for(i=1;i<=NF;i+=4) a+=$i; for(i=2;i<=NF;i+=4) c+=$i; for(i=3;i<=NF;i+=4) g+=$i; for(i=4;i<=NF;i+=4) t+=$i; print a\"\\t\"c\"\\t\"g\"\\t\"t}'",
    shQuote(c_file)
  )
  counts <- fread(cmd = awk_cmd, col.names = c("A", "C", "G", "T"))

  if (nrow(counts) != nrow(maf)) {
    stop(sprintf("Row count mismatch between %s (%d) and %s (%d)", m_file, nrow(maf), c_file, nrow(counts)))
  }

  dt <- cbind(maf, counts)
  
  # Restrict to called SNPs only
  dt <- dt[snps, on = .(chromo, position), nomatch = 0]
  dt[, `:=`(pop = target_pop, era = target_era, pop_era = pop_era)]
  
  # 3. Calculate Reference Bias (Ref Count / (Ref Count + Alt Count))
  # Identify the non-reference allele between major and minor
  dt[, alt_allele := ifelse(ref == major, minor, major)]
  dt[, ref_count := ifelse(ref == "A", A, ifelse(ref == "C", C, ifelse(ref == "G", G, ifelse(ref == "T", T, 0))))]
  dt[, alt_count := ifelse(alt_allele == "A", A, ifelse(alt_allele == "C", C, ifelse(alt_allele == "G", G, ifelse(alt_allele == "T", T, 0))))]
  dt[, total_depth := ref_count + alt_count]
  
  # Filter for sites with valid depth to compute allele balance
  dt_valid <- dt[total_depth > 0]
  dt_valid[, ref_ratio := ref_count / total_depth]
  
  mean_ref_ratio <- mean(dt_valid$ref_ratio, na.rm = TRUE)
  
  # ---------------------------------------------------------------------------
  # 4. Calculate Strand Bias from snpStat file (with Tarone's Criterion)
  # ---------------------------------------------------------------------------
  # ANGSD's snpStat.gz reports forward (+)/reverse (-) read counts for Major
  # and Minor alleles alongside SB1:SB2:SB3, where SB3 is the two-sided
  # Fisher's exact test p-value on the 2x2 strand-by-allele contingency table.
  #
  # Why we apply Tarone's criterion (Tarone 1990, Biometrics):
  #   1. Discrete test power: Because read counts are integers, SNPs with low
  #      minor-allele or single-strand coverage cannot reach p < 0.05 even under
  #      100% strand asymmetry.
  #   2. Monomorphic sites: ANGSD outputs SB3 = 1.000000 for global SNPs that
  #      are monomorphic within a given cohort (+Minor + -Minor == 0).
  #   Including these untestable/monomorphic sites in the denominator artificially
  #   deflates the % significant strand-bias rate in lower-depth cohorts (e.g.,
  #   historic samples). Conditioning on observed marginal totals to compute the
  #   minimum achievable p-value (p_min) and keeping only loci with p_min < 0.05
  #   removes both monomorphic and underpowered sites without biasing the test.
  # ---------------------------------------------------------------------------
  s_match <- snpstat_files[basename(snpstat_files) == paste0(pop_era, ".snpStat.gz")]
  if (length(s_match) != 1 || !file.exists(s_match)) {
    stop(sprintf("Could not uniquely match snpStat file for %s", pop_era))
  }
  s_file <- s_match[1]

  # ANGSD -doSNPstat 1 separates "+Major +Minor -Major -Minor" and their four
  # integer counts (e.g., "15 0 20 1") with spaces rather than tabs on every line,
  # while separating all other fields with tabs. Rebuilding every record with
  # OFS="\t" ({$1=$1}) converts all mixed space/tab delimiters into 11 uniform
  # tab-delimited columns so fread() can select and validate columns by name.
  snpstat_cmd <- sprintf("zcat %s | awk 'BEGIN{OFS=\"\\t\"} {$1=$1} 1'", shQuote(s_file))
  snpstat <- fread(
    cmd = snpstat_cmd,
    select = c("Chromo", "Position", "+Major", "+Minor", "-Major", "-Minor", "SB1:SB2:SB3")
  )
  setnames(
    snpstat,
    c("Chromo", "Position", "+Major", "+Minor", "-Major", "-Minor", "SB1:SB2:SB3"),
    c("chromo", "position", "plus_maj", "plus_min", "minus_maj", "minus_min", "sb_info")
  )

  # Subset to the global set of called SNPs before parsing SB1:SB2:SB3
  dt_sb <- snpstat[snps, on = .(chromo, position), nomatch = 0]
  sb3_raw <- tstrsplit(dt_sb$sb_info, ":", fixed = TRUE)[[3]]
  dt_sb[, sb_pval := suppressWarnings(as.numeric(sb3_raw))]

  # Marginal totals of each site's 2x2 table:
  #   R1 = total forward (+) reads, R2 = total reverse (-) reads, C2 = total minor reads
  dt_sb[, `:=`(
    R1 = plus_maj + plus_min,
    R2 = minus_maj + minus_min,
    C2 = plus_min + minus_min
  )]

  # Under fixed margins, forward minor reads follow Hypergeom(R1, R2, C2) on [b_min, b_max].
  # Because the hypergeometric PMF is unimodal, the least probable table (and thus p_min)
  # always lies at one of the two boundary counts b_min or b_max.
  dt_sb[, `:=`(
    b_min = pmax(0L, C2 - R2),
    b_max = pmin(R1, C2)
  )]
  dt_sb[, `:=`(
    p1 = dhyper(b_min, m = R1, n = R2, k = C2),
    p2 = dhyper(b_max, m = R1, n = R2, k = C2)
  )]

  # Edge cases for p_min:
  #   - b_min == b_max (monomorphic or single-strand sites): only 1 table exists -> p_min = 1.0
  #   - p1 == p2 (symmetric margins): both tails contribute to two-sided p -> p_min = p1 + p2
  dt_sb[, p_min := fifelse(
    b_min == b_max, 1.0,
    fifelse(abs(p1 - p2) <= 1e-12 * pmax(p1, p2), p1 + p2, pmin(p1, p2))
  )]

  # Retain only finite, Tarone-testable loci (p_min < 0.05) for the denominator
  dt_sb <- dt_sb[is.finite(sb_pval) & is.finite(p_min) & p_min < 0.05]
  pct_sb_significant <- if (nrow(dt_sb) > 0) mean(dt_sb$sb_pval < 0.05) * 100 else 0

  results_list[[i]] <- data.table(
    pop = target_pop,
    era = target_era,
    pop_era = pop_era,
    mean_ref_ratio = mean_ref_ratio,
    n_sb_testable = nrow(dt_sb),
    pct_strand_bias = pct_sb_significant
  )
  
  maf_dt_list[[i]] <- dt_valid[, .(pop, era, pop_era, ref_ratio)]
}

summary_dt <- rbindlist(results_list)
all_sites_dt <- rbindlist(maf_dt_list)

# Write summary output table
fwrite(summary_dt, alist$out_tsv, sep = "\t")

# Generate 2-Panel QC Bias Plot (Reference Bias + Strand Bias)
p_ref <- ggplot(all_sites_dt, aes(x = ref_ratio, fill = era)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  theme_minimal() +
  facet_wrap(~pop) +
  labs(title = "Reference Bias (Allele Balance)", x = "Ref / (Ref + Alt) Ratio", y = "Density", fill = "Era")

# Bar plot of % significant strand-bias loci among Tarone-testable SNPs (p_min < 0.05)
p_sb <- ggplot(summary_dt, aes(x = pop, y = pct_strand_bias, fill = era)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = sprintf("%.2f%%", pct_strand_bias)),
    position = position_dodge(width = 0.9),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, NA)) +
  theme_minimal() +
  labs(
    title = "Strand Bias (% Testable Loci with Fisher's Exact p < 0.05)",
    subtitle = "Denominator restricted to Tarone-testable SNPs (p_min < 0.05; excludes monomorphic/low-power sites)",
    x = "Region/Population",
    y = "% Significant Testable Loci",
    fill = "Era"
  )

png(alist$out_png, width = 10, height = 8, units = "in", res = 300, type = "cairo")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p_ref, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sb,  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()


# Individual QC: Mean Depth vs. Observed Heterozygosity (Ho)
ind_depth_dt <- rbindlist(ind_depth_list)

# Parse 1D SFS (<sample>.est.ml: [hom_anc, het, hom_der])
sfs_dt <- rbindlist(lapply(ind_sfs_files, function(f) {
  s_name <- sub("\\.est\\.ml$", "", basename(f))
  vals <- scan(f, quiet = TRUE)
  ho <- if (length(vals) >= 3 && sum(vals) > 0) vals[2] / sum(vals) else NA_real_
  data.table(sample = s_name, Ho = ho)
}))

ind_qc_dt <- merge(ind_depth_dt, sfs_dt, by = "sample", all.x = TRUE)
fwrite(ind_qc_dt, alist$out_ind_tsv, sep = "\t")

p_het <- ggplot(ind_qc_dt, aes(x = mean_depth, y = Ho, color = era, shape = pop)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(aes(group = era, fill = era), method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.7, alpha = 0.15) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Individual Observed Heterozygosity (Ho) vs. Mean Callable Depth",
    x = "Mean Depth Across Callable Sites (X)",
    y = "Observed Heterozygosity (Ho)",
    color = "Era",
    fill = "Era",
    shape = "Region/Population"
  )

ggsave(alist$out_ind_png, plot = p_het, width = 7.5, height = 5, dpi = 300)