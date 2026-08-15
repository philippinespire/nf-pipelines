process ANGSD_FST {
    tag "${pop}_hist_vs_mod"
    
    publishDir "${params.outdir}/fst",                 mode: 'copy', pattern: "*.{txt,2dsfs}"
    publishDir "${params.outdir}/large_data/fst",      mode: 'copy', pattern: "*.fst.*"

    input:
    tuple val(pop), path(hist_saf), path(hist_idx), path(hist_pos), path(mod_saf), path(mod_idx), path(mod_pos)
    tuple path(sites_pos), path(sites_bin), path(sites_idx)

    output:
    path "${pop}_hist_vs_mod.2dsfs",          emit: sfs_2d
    path "${pop}_hist_vs_mod.global_fst.txt", emit: global_fst
    path "${pop}_hist_vs_mod.fst.gz",         emit: fst_gz
    path "${pop}_hist_vs_mod.fst.idx",        emit: fst_idx

    script:
    """
    # 1. Estimate 2D-SFS using only LD-pruned neutral sites
    realSFS ${hist_idx} ${mod_idx} \
        -sites ${sites_pos} \
        -P ${task.cpus} > ${pop}_hist_vs_mod.2dsfs

    # 2. Index FST per-site
    realSFS fst index ${hist_idx} ${mod_idx} \
        -sfs ${pop}_hist_vs_mod.2dsfs \
        -sites ${sites_pos} \
        -out ${pop}_hist_vs_mod \
        -P ${task.cpus}

    # 3. Calculate global weighted and unweighted FST statistics
    realSFS fst stats ${pop}_hist_vs_mod.fst.idx > ${pop}_hist_vs_mod.global_fst.txt
    """
}