process ANGSD_DIVERSITY {
    tag "${pop}_${era}"
    publishDir "${params.outdir}/diversity",          mode: 'copy', pattern: "*.{sfs,pestPG}"
    publishDir "${params.outdir}/large_data/thetas", mode: 'copy', pattern: "*.thetas*"
    
    input:
    tuple val(pop), val(era), path(saf), path(saf_idx), path(saf_pos)
    path neutral_sites  // [pos, bin, idx] bundle from PREPARE_DIVERSITY_SITES

    output:
    path "${pop}_${era}.sfs", emit: sfs
    path "${pop}_${era}.pestPG", emit: pestPG
    path "${pop}_${era}.thetas.idx", emit: thetas_idx
    path "${pop}_${era}.thetas.gz", emit: thetas

    script:
    """
    # 1. Generate Unpruned SFS using neutral callable sites (includes monomorphic loci)
    realSFS ${pop}_${era}.saf.idx -sites ${neutral_sites[0]} -P ${task.cpus ?: 32} -fold 1 > ${pop}_${era}.sfs

    # 2. Calculate Pi and Theta using the Neutral SFS
    realSFS saf2theta ${pop}_${era}.saf.idx -sites ${neutral_sites[0]} -sfs ${pop}_${era}.sfs -fold 1 -P ${task.cpus ?: 32} -outname ${pop}_${era}
    thetaStat do_stat ${pop}_${era}.thetas.idx -outnames ${pop}_${era}
    """
}