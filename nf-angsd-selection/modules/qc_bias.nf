process EVALUATE_QC_BIAS {
    tag "Evaluate Historic vs Modern Biases"
    publishDir "${params.outdir}/qc_bias", mode: 'copy'

    input:
    path maf_files    // Collected .mafs.gz files from ANGSD_GL_POP
    path count_files  // Collected .counts.gz files from ANGSD_GL_POP
    path snpstat_files  // Collected .snpStat.gz files
    path snps_pos     // PREPARE_SITE_SETS.out.all_snps
    path samplesheet

    output:
    path "qc_bias_summary.tsv", emit: summary
    path "qc_bias_plots.png"  , emit: plot

    script:
    """
    Rscript ${projectDir}/scripts/compare_qc_bias.R \
        --mafs="${maf_files.join(',')}" \
        --counts="${count_files.join(',')}" \
        --snpstat="${snpstat_files.join(',')}" \
        --snps="${snps_pos}" \
        --samplesheet="${samplesheet}" \
        --out_tsv="qc_bias_summary.tsv" \
        --out_png="qc_bias_plots.png"
    """
}