process CALC_IND_HET {
    tag "${sample}"

    input:
    tuple val(sample), val(era), val(region), path(bam)
    tuple path(ref), path(ref_fai)
    tuple path(callable_pos), path(callable_bin), path(callable_idx)

    output:
    path "${sample}.est.ml", emit: sfs

    script:
    """
    # Ensure BAM index is visible if not auto-staged
    [ -e "${bam}.bai" ] || ln -s "\$(readlink -f ${bam}).bai" "${bam}.bai" 2>/dev/null || true

    angsd -i ${bam} \
        -ref "${ref}" \
        -anc "${ref}" \
        -dosaf 1 \
        -GL 1 \
        -minMapQ 25 -minQ 30 \
        -uniqueOnly 1 -remove_bads 1 \
        -sites ${callable_pos} \
        -noTrans 1 \
        -P ${task.cpus ?: 2} \
        -out ${sample}

    realSFS ${sample}.saf.idx -P ${task.cpus ?: 2} > ${sample}.est.ml
    """
}

process EVALUATE_QC_BIAS {
    tag "Evaluate Historic vs Modern Biases"
    publishDir "${params.outdir}/qc_bias", mode: 'copy'

    input:
    path maf_files      // Collected .mafs.gz files from ANGSD_GL_POP
    path count_files    // Collected .counts.gz files from ANGSD_GL_POP
    path snpstat_files  // Collected .snpStat.gz files
    path ind_sfs_files  // Collected per-individual .est.ml files from CALC_IND_HET
    path snps_pos       // PREPARE_SITE_SETS.out.all_snps
    path samplesheet

    output:
    path "qc_bias_summary.tsv", emit: summary
    path "ind_het_depth.tsv"   , emit: ind_het
    path "qc_bias_plots.png"  , emit: plot
    path "ind_het_vs_depth.png", emit: ind_plot

    script:
    """
    Rscript ${projectDir}/scripts/compare_qc_bias.R \
        --mafs="${maf_files.join(',')}" \
        --counts="${count_files.join(',')}" \
        --snpstat="${snpstat_files.join(',')}" \
        --ind_sfs="${ind_sfs_files.join(',')}" \
        --snps="${snps_pos}" \
        --samplesheet="${samplesheet}" \
        --out_tsv="qc_bias_summary.tsv" \
        --out_ind_tsv="ind_het_depth.tsv" \
        --out_png="qc_bias_plots.png" \
        --out_ind_png="ind_het_vs_depth.png"
    """
}
