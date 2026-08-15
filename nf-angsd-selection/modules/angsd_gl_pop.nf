// Identify loci passing minind_frac and maxdepth_mult per population/era subset
process QC_POP_SITES {
    tag "${pop}_${era}"

    input:
    tuple val(pop), val(era), path(bamlist)
    path bed_sites // [snps, bin, idx] from INDEX_BED_SITES

    output:
    path "${pop}_${era}.passing.pos", emit: passing_pos

    script:
    """
    # 1. Total target length from BED file
    target_bp=\$(awk 'BEGIN{s=0} {s += \$3 - \$2} END{print s}' ${params.bed_file})

    # 2. Calculate subset-specific total mapped bases (historical vs modern)
    total_bases=0
    while read -r bam; do
        n_reads=\$(samtools idxstats "\$bam" | awk '{s+=\$3} END{print s}')
        read_len=\$(samtools view "\$bam" | head -n 1000 | awk 'NR>0 {s+=length(\$10); c++} END{if(c>0) print int(s/c); else print 100}')
        total_bases=\$(awk -v tb="\$total_bases" -v nr="\$n_reads" -v rl="\$read_len" 'BEGIN{print tb + (nr * rl)}')
    done < ${bamlist}

    # 3. Calculate subset-specific minInd and max_depth
    n_ind=\$(wc -l < ${bamlist})
    min_ind=\$(awk -v n="\$n_ind" -v f="${params.min_ind_frac}" 'BEGIN{print int(n * f + 0.999)}')
    
    mean_combined_depth=\$(awk -v tb="\$total_bases" -v tbp="\$target_bp" 'BEGIN{print tb / tbp}')
    max_depth=\$(awk -v md="\$mean_combined_depth" -v m="${params.max_depth_mult}" 'BEGIN{print int(md * m + 0.999)}')

    # 4. Run ANGSD QC pass for this specific population/era
    angsd -bam ${bamlist} \
        -GL 1 -doMajorMinor 4 -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -minInd \$min_ind \
        -setmaxdepth \$max_depth \
        -sites ${bed_sites[0]} \
        -out ${pop}_${era}_qc \
        -P ${task.cpus}

    # 5. Extract loci that passed these subset-specific thresholds
    zcat ${pop}_${era}_qc.mafs.gz | awk 'NR>1 {print \$1"\t"\$2}' > ${pop}_${era}.passing.pos
    """
}

// Intersect passing sites across both historical and modern population/era subsets (Venn Diagram Overlap)
process INTERSECT_POP_SITES {
    tag "Loci Intersection Across Populations"
    publishDir "${params.outdir}/sites", mode: 'copy'

    input:
    path pos_files // Collect of all passing.pos files

    output:
    path "pop_intersect.pos"    , emit: pos
    path "pop_intersect.pos.bin", emit: bin
    path "pop_intersect.pos.idx", emit: idx

    script:
    """
    num_files=\$(ls -1 ${pos_files} | wc -l)

    # Retain ONLY loci that pass QC filters in EVERY historical and modern subset
    cat ${pos_files} | sort | uniq -c | awk -v n="\$num_files" '\$1 == n {print \$2"\t"\$3}' > pop_intersect.pos

    angsd sites index pop_intersect.pos
    """
}

// Final Population Execution using the shared intersected site index
// run ANGSD on each population and era combination
process ANGSD_GL_POP {
    tag { "${pop}_${era}" }
    publishDir "${params.outdir}/large_data/saf",  mode: 'copy', pattern: "*.saf*"
    publishDir "${params.outdir}/large_data/mafs", mode: 'copy', pattern: "*.mafs.gz"

    input:
    tuple val(pop), val (era), path(bamlist)
    path intersect_snps // pop_intersect.pos
    path intersect_bin  // pop_intersect.pos.bin
    path intersect_idx  // pop_intersect.pos.idx
    path regions

    output:
    tuple val(pop), val(era), path("${pop}_${era}.saf.gz"), path("${pop}_${era}.saf.idx"), path("${pop}_${era}.saf.pos.gz"), emit: saf_files
    path "${pop}_${era}.mafs.gz", emit: mafs

    script:
    """
        angsd -bam ${bamlist} \
        -doSaf 1 \
        -GL 1 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 \
        -P ${task.cpus} \
        -ref ${params.reference} \
        -anc ${params.reference} \
        -sites ${intersect_snps} \
        -rf ${regions} \
        -out ${pop}_${era} \
        -noTrans 1
    """
}
