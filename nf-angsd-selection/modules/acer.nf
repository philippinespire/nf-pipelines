process RUN_ACER {
    tag "ACER Selection Scan"
    publishDir "${params.outdir}/selection", mode: 'copy'

    input:
    path mafs           // All ${region}_${era}.mafs.gz files from ANGSD_GL_POP
    path var_sites      // Polymorphic sites index (sites.snps)
    path helpers        // acer_helpers.R
    path run_script     // run_acer.R
    path samplesheet    // To map historic/modern pairs per region

    output:
    path "iteration_summary.tsv"     , emit: summary
    path "final_test_results.tsv"    , emit: chisq_results
    path "ne_bootstrap.tsv"          , emit: ne_boot
    path "*.png"                     , optional: true, emit: plots

    script:
    """
    # 1. Intersect population MAF files with polymorphic SNP list to discard monomorphic sites
    for maf in *.mafs.gz; do
        zcat \$maf | awk 'NR==FNR {a[\$1"\t"\$2]; next} FNR==1 || (\$1"\t"\$2) in a' ${var_sites} - | gzip > var_\${maf}
    done

    # 2. Parse sample CSV to construct --hist_mafs, --mod_mafs, and --region_names CLI arguments
    read -r HIST_MAFS MOD_MAFS REGIONS <<< \$(Rscript -e '
        df <- read.csv("${samplesheet}", stringsAsFactors = FALSE)
        df\$region <- trimws(df\$region)
        df\$era_norm <- tolower(trimws(df\$era))
        hist_set <- c("historic", "historical", "hist", "era1")
        mod_set  <- c("modern", "mod", "era2")

        regs <- unique(df\$region)
        paired_regs <- Filter(function(r) {
            any(df\$region == r & df\$era_norm %in% hist_set) &&
            any(df\$region == r & df\$era_norm %in% mod_set)
        }, regs)

        if (length(paired_regs) == 0) {
            stop("No regions with both historic and modern samples found in samplesheet.")
        }

        hist_files <- sapply(paired_regs, function(r) {
            era <- df\$era_norm[df\$region == r & df\$era_norm %in% hist_set][1]
            paste0("var_", r, "_", era, ".mafs.gz")
        })
        mod_files <- sapply(paired_regs, function(r) {
            era <- df\$era_norm[df\$region == r & df\$era_norm %in% mod_set][1]
            paste0("var_", r, "_", era, ".mafs.gz")
        })

        cat(paste(hist_files, collapse=","), paste(mod_files, collapse=","), paste(paired_regs, collapse=","), "\n")
    ')

    # 3. Execute ACER iterative selection scan
    Rscript ${run_script} \\
        --hist_mafs="\$HIST_MAFS" \\
        --mod_mafs="\$MOD_MAFS" \\
        --region_names="\$REGIONS" \\
        --out_dir="." \\
        --helpers="${helpers}" \\
        --generations=${params.generations} \\
        --fdr_cutoff=${params.fdr_cutoff} \\
        --max_rounds=${params.max_rounds} \\
        --n_boot=${params.n_boot} \\
        --min_ind=${params.min_ind}
    """
}
