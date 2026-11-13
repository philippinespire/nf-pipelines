#!/usr/bin/env Rscript

# goal is to calculate Reference Bias (Ref / (Ref + Alt)), and Strand Bias (p < 0.05 fraction) across eras.

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

results_list <- list()
maf_dt_list  <- list()

# Read the target SNP list
snps <- fread(alist$snps, select = 1:2, col.names = c("chromo", "position"))

for (i in seq_along(maf_files)) {
  m_file <- maf_files[i]
  
  # Extract population and era from filename (e.g., Bali_historic.mafs.gz)
  base_name <- basename(m_file)
  pop_era   <- gsub("\\.mafs\\.gz$", "", base_name)
  era       <- ifelse(grepl("_(historic|historical|hist|era1)$", pop_era, ignore.case = TRUE), "Historic", "Modern")
  pop       <- sub("_(historic|historical|hist|era1|modern|mod|era2)$", "", pop_era, ignore.case = TRUE)
  
  # Match corresponding counts file by exact basename (Nextflow .collect() order is not guaranteed)
  c_match <- count_files[basename(count_files) == paste0(pop_era, ".counts.gz")]
  if (length(c_match) != 1 || !file.exists(c_match)) {
    stop(sprintf("Could not uniquely match counts file for %s", pop_era))
  }
  c_file <- c_match[1]

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
  dt[, `:=`(pop = pop, era = era, pop_era = pop_era)]
  
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
  
  # 4. Calculate Strand Bias from snpStat file
  s_match <- snpstat_files[basename(snpstat_files) == paste0(pop_era, ".snpStat.gz")]
  if (length(s_match) != 1 || !file.exists(s_match)) {
    stop(sprintf("Could not uniquely match snpStat file for %s", pop_era))
  }
  s_file <- s_match[1]

  snpstat <- fread(cmd = paste("zcat", shQuote(s_file)), select = c("Chromo", "Position", "SB1:SB2:SB3"))
  setnames(snpstat, c("Chromo", "Position", "SB1:SB2:SB3"), c("chromo", "position", "sb_info"))

  # Subset to called SNPs before parsing SB1:SB2:SB3
  dt_sb <- snpstat[snps, on = .(chromo, position), nomatch = 0]
  sb3_raw <- tstrsplit(dt_sb$sb_info, ":", fixed = TRUE)[[3]]
  dt_sb[, sb_pval := suppressWarnings(as.numeric(sb3_raw))]
  dt_sb <- dt_sb[is.finite(sb_pval)]
  pct_sb_significant <- if (nrow(dt_sb) > 0) mean(dt_sb$sb_pval < 0.05) * 100 else 0
  
  results_list[[i]] <- data.table(
    pop = pop,
    era = era,
    pop_era = pop_era,
    mean_ref_ratio = mean_ref_ratio,
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
  labs(title = "Strand Bias (% Loci p < 0.05)", x = "Region/Population", y = "% Significant Loci", fill = "Era")

png(alist$out_png, width = 10, height = 8, units = "in", res = 300, type = "cairo")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(p_ref, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sb,  vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()