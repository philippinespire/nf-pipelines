// 1. Prepare Neutral Callable Sites for Diversity (Includes Monomorphic Sites)
process PREPARE_DIVERSITY_SITES {
    tag "Neutral Callable Sites (Diversity)"
    publishDir "${params.outdir}/filtered_sites", mode: 'copy'

    input:
    path callable_bed  // params.bed_file
    path chisq_results   // ACER output or []

    output:
    path "diversity_neutral.pos"    , emit: pos
    path "diversity_neutral.pos.bin", emit: bin
    path "diversity_neutral.pos.idx", emit: idx

    script:
    """
    # Convert BED file (chr start end) to 2-column POS file (chr pos)
    awk 'BEGIN{OFS="\t"} {for(i=\$2+1; i<=\$3; i++) print \$1, i}' ${callable_bed} > all_callable.pos

    # Exclude candidate selection loci if present
    if [ -f "${chisq_results}" ] && [ -s "${chisq_results}" ]; then
        awk -v fdr="${params.fdr_cutoff}" 'NR>1 && \$NF <= fdr {print \$1"\t"\$2}' ${chisq_results} > selected_coords.txt
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_coords.txt all_callable.pos > diversity_neutral.pos
    else
        cp all_callable.pos diversity_neutral.pos
    fi

    # Index for ANGSD
    angsd sites index diversity_neutral.pos
    """
}

// 2. Prepare Neutral & LD-Pruned Sites for PCA & Admixture (Polymorphic Only)
process PREPARE_PCA_SITES {
    tag "Neutral & Pruned Loci (PCA/Admix)"
    publishDir "${params.outdir}/filtered_sites", mode: 'copy'

    input:
    path poly_snps   // sites.snps (polymorphic sites)
    path pruned_pos  // LD-pruned sites or []
    path chisq_results // ACER output or []

    output:
    path "pca_neutral.pos", emit: pos

    script:
    """
    # Base set: Use LD-pruned sites if available, otherwise use all polymorphic sites
    if [ -f "${pruned_pos}" ] && [ -s "${pruned_pos}" ]; then
        cp ${pruned_pos} base_pca.pos
    else
        cp ${poly_snps} base_pca.pos
    fi

    # Exclude candidate selection loci if present
    if [ -f "${chisq_results}" ] && [ -s "${chisq_results}" ]; then
        awk -v fdr="${params.fdr_cutoff}" 'NR>1 && \$NF <= fdr {print \$1"\t"\$2}' ${chisq_results} > selected_coords.txt
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_coords.txt base_pca.pos > pca_neutral.pos
    else
        cp base_pca.pos pca_neutral.pos
    fi
    """
}

process PREPARE_SITE_SETS {
    tag "Generate Complete Site Set Library"
    publishDir "${params.outdir}/sites", mode: 'copy'

    input:
    path pop_intersect  // pop_intersect.pos (Venn diagram output)
    path poly_snps     // sites.snps from ANGSD_GL_ALL
    path pruned_pos    // merged_sites from LD_PRUNE (or [])
    path chisq_results // ACER output (or [])

    output:
    path "all_callable.pos"             , emit: all_callable
    path "all_snps.pos"                 , emit: all_snps
    path "selected_loci.pos"            , emit: selected
    path "callable_neutral.pos"         , emit: callable_neutral_pos
    path "callable_neutral.pos.bin"     , emit: callable_neutral_bin
    path "callable_neutral.pos.idx"     , emit: callable_neutral_idx
    path "snps_neutral.pos"             , emit: snps_neutral
    path "snps_pruned.pos"              , emit: snps_pruned
    path "snps_pruned_neutral.pos"      , emit: snps_pruned_neutral_pos
    path "snps_pruned_neutral.pos.bin"  , emit: snps_pruned_neutral_bin
    path "snps_pruned_neutral.pos.idx"  , emit: snps_pruned_neutral_idx

    script:
    """
    module load container_env angsd

    cp ${pop_intersect} all_callable.pos
    awk 'NR==FNR {pop[\$1"\t"\$2]; next} (\$1"\t"\$2) in pop' all_callable.pos ${poly_snps} > all_snps.pos

    if [ -f "${pruned_pos}" ] && [ -s "${pruned_pos}" ]; then
        awk 'NR==FNR {pop[\$1"\t"\$2]; next} (\$1"\t"\$2) in pop' all_callable.pos ${pruned_pos} > snps_pruned.pos
    else
        cp all_snps.pos snps_pruned.pos
    fi

    if [ -f "${chisq_results}" ] && [ -s "${chisq_results}" ]; then
        awk -v fdr="${params.fdr_cutoff}" 'NR>1 && \$NF <= fdr {print \$1"\t"\$2}' ${chisq_results} > selected_loci.pos
    else
        touch selected_loci.pos
    fi

    if [ -s selected_loci.pos ]; then
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos all_callable.pos > callable_neutral.pos
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos all_snps.pos     > snps_neutral.pos
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos snps_pruned.pos   > snps_pruned_neutral.pos
    else
        cp all_callable.pos callable_neutral.pos
        cp all_snps.pos     snps_neutral.pos
        cp snps_pruned.pos  snps_pruned_neutral.pos
    fi

    # Index position lists for realSFS site filtering
    angsd sites index callable_neutral.pos
    angsd sites index snps_pruned_neutral.pos
    """
}