process PCANGSD {
    // calculates pca covariance matrix and admixture proportions
    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path (all_beagle)

    output:
    path "${params.species}.cov", emit: pcangsd_cov
	path "${params.species}.*.Q", emit: q_matrix

    script:
    """
    # if <100 sites, run admixture in one batch to avoid errors (5 or 10 is default to speed up computation)
    n_sites=\$(( \$(gzip -dc ${all_beagle} | wc -l) - 1 ))
    admix_batch_opt=""
    if [ "\${n_sites}" -lt 100 ]; then
        admix_batch_opt="--admix-batch 1"
    fi

    pcangsd -b ${all_beagle} --admix \${admix_batch_opt} -t ${task.cpus ?: 8} -o ${params.species}
    """
}

process PLOT_PCANGSD {

    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path (pcangsd_cov)
    path samplesheet_file

    output:
    path "${params.species}.pcangsd.plot.png"

    script:
    """
    Rscript ${projectDir}/scripts/pcangsd.R \\
        --cov ${pcangsd_cov} \\
        --samplesheet ${samplesheet_file} \\
        --species ${params.species} \\
        --out ${params.species}.pcangsd.plot.png
    """
}

process PLOT_ADMIXTURE {
    
    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path q_file
    path metadata

    output:
    path "${params.species}.admixture.png"

    script:
    """
    # Runs the provided R script using the matrix and metadata
    Rscript ${projectDir}/scripts/plot_admixture.R \
        ${q_file} \
        ${metadata} \
        ${params.species}.admixture.png
    """
}