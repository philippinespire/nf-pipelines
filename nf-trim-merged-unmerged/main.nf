nextflow.enable.dsl=2

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt"
params.indir        = "${projectDir}/data/symlinks"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/<reference>.fasta"
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed"
params.reference_prefix         = params.reference.tokenize('/').last().replaceAll(/\.(fa|fasta|fna)$/, '')
params.historical_era           = "historical"
params.run_historical_fastqc    = false
params.run_historical_mapdamage = false
params.use_historical_rescaled  = false
params.run_historical_amber     = false
params.run_modern_amber         = false
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/archive/carpenterlab/pire/softwares/AMBERv2/AMBER" 
params.bwa_threads  = 4
params.bam_q        = 1 // Mapping quality. Currently set to 1 simply to remove unmapped reads. 
params.trimlength   = 85 // the default trim length if no historical reads are available to calculate it from

def resolve_reads = { String sample_id ->
    def r1 = file("${params.indir}/${sample_id}_R1.fastq.gz")
    def r2 = file("${params.indir}/${sample_id}_R2.fastq.gz")

    if( !r1.exists() ) error "R1 file not found for ${sample_id}: ${r1}"
    if( !r2.exists() ) error "R2 file not found for ${sample_id}: ${r2}"

    return [r1, r2]
}

def samplesheet_file = file(params.samplesheet)
def legacy_samples_file = file(params.samples_file)

// --- Input Channels ---
// Sample metadata is taken from a CSV when available.
// Fallback: legacy one-column sample list, which defaults all samples to modern.
sample_metadata_ch = samplesheet_file.exists() ?
    Channel
        .fromPath(samplesheet_file, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample.trim(), row.era?.trim()?.toLowerCase() ?: 'modern') }
        .filter { sample_id, era -> sample_id }
    :
    Channel
        .fromPath(legacy_samples_file, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }
        .map { sample_id -> tuple(sample_id, 'modern') }

modern_samples_ch = sample_metadata_ch
    .filter { sample_id, era -> era != params.historical_era }
    .map { sample_id, era -> tuple(sample_id, resolve_reads(sample_id)) }

historical_samples_ch = sample_metadata_ch
    .filter { sample_id, era -> era == params.historical_era }
    .map { sample_id, era ->
        def reads = resolve_reads(sample_id)
        tuple(sample_id, reads[0], reads[1])
    }

// For historical reads: extract paths for trim length calculation
    historical_read_paths_ch = historical_samples_ch
        .flatMap { sample_id, r1, r2 -> [r1.toString(), r2.toString()] }
        .collect()

// --- Processes ---


process FASTP_MERGED {
    // takes paired-end raw sequencing reads, cleans them up (quality filtering, adapter trimming, poly-G tail removal), 
    // and attempts to merge overlapping forward and reverse reads into single, longer reads.
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.fastq.gz"), emit: merged
    tuple val(sample_id), path("${sample_id}_R1_unmerged.fq.gz"), path("${sample_id}_R2_unmerged.fq.gz"), emit: unmerged
    path "*.json"
    path "*.html"

    script:
    """
    fastp \
      -i ${reads[0]} -I ${reads[1]} \
      -p -c --trim_poly_g --merge \
      --merged_out=${sample_id}_trimmed_merged.fastq.gz \
      -o ${sample_id}_R1_unmerged.fq.gz -O ${sample_id}_R2_unmerged.fq.gz \
      -h ${sample_id}_fastp_report.html -j ${sample_id}_fastp_report.json \
      -R "${sample_id}" -w ${task.cpus} -l 30 --overlap_diff_limit 1 --overlap_len_require 11
    """
}

process SPLIT_MERGED {
    // processes a FASTQ file and splits any read that is >trim_len right down the middle, 
    // turning a single long read into two "pseudo-paired" reads. 
    // Any read that is already shorter than your target length is left completely untouched.
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq), val(trim_len)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${trim_len}.fastq.gz")

    script:
    """
    echo "Looking for script at: ${params.split_script}"
    ls -l ${params.split_script}
    bash ${params.split_script} ${merged_fq} ${sample_id}_trimmed_merged.L${trim_len}.fastq.gz ${trim_len}
    """
}

process QC_MERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -o . -t ${task.cpus} --extract ${merged_fq}
    """
}

process FASTP_UNMERGED {
    // filters, trims, and generates quality reports for the unmerged paired-end reads, keeping them separate
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fq.gz"), path("${sample_id}_R2.trimmed.fq.gz")
    path "*.json"
    path "*.html"

    script:
    """
    fastp \
      -i ${r1} -I ${r2} \
      -p -c --trim_poly_g \
      -o ${sample_id}_R1.trimmed.fq.gz -O ${sample_id}_R2.trimmed.fq.gz \
      -h ${sample_id}_unmerged_fastp_report.html \
      -j ${sample_id}_unmerged_fastp_report.json \
      -R "${sample_id}" -w ${task.cpus} -l 30
    """
}

process SEQTK_TRIM {
    // takes unmerged paired-end reads (the ones that didn't overlap enough to be combined) 
    // and cleanly truncates them down to a strict maximum length (trim_len) using seqtk
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), val(trim_len)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.L${trim_len}.fastq.gz"), path("${sample_id}_R2_trimmed.L${trim_len}.fastq.gz")

    script:
    """
    seqtk trimfq -L ${trim_len} ${r1} | gzip > ${sample_id}_R1_trimmed.L${trim_len}.fastq.gz
    seqtk trimfq -L ${trim_len} ${r2} | gzip > ${sample_id}_R2_trimmed.L${trim_len}.fastq.gz
    """
}

process QC_UNMERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -o . -t ${task.cpus} --extract ${r1}
    fastqc -o . -t ${task.cpus} --extract ${r2}
    """
}

process PREP_REFERENCE_REPEAT {
    publishDir "${params.outdir}/data/reference", mode: 'copy'

    input:
    path ref
    path repeat_bed

    output:
    path "${params.reference_prefix}.repma.angsd.txt", emit: sites
    path "${params.reference_prefix}.repma.angsd.txt.idx", emit: sites_idx
    path "${params.reference_prefix}.repma.angsd.txt.bin", emit: sites_bin
    path "${params.reference_prefix}.regions", emit: regions
    path "${params.reference_prefix}.chrs", emit: chrs

    script:
    """
    awk '{print \$1"\t"\$2+1"\t"\$3}' ${repeat_bed} > ${params.reference_prefix}.repma.angsd.txt
    angsd sites index ${params.reference_prefix}.repma.angsd.txt
    cut -f1 ${params.reference_prefix}.repma.angsd.txt | awk '!seen[\$0]++' | awk '{print \$0 ":"}' > ${params.reference_prefix}.regions
    cut -f1 ${params.reference_prefix}.repma.angsd.txt | sort | uniq > ${params.reference_prefix}.chrs
    """
}

process CALC_HISTORICAL_TRIMLEN {
    publishDir "${params.outdir}/data/stats", mode: 'copy'

    input:
    val read_paths

    output:
    stdout emit: trim_len
    path "historical_trimlength.txt", emit: trim_file

    script:
    def reads_list = read_paths.collect { "'${it}'" }.join(', ')
    """
    python3 - <<'PY'
    import gzip

    # List of gzipped FASTQ file paths dynamically injected by Nextflow
    files = [${reads_list}]
    
    # Counters to track the total number of nucleotides (bases) and individual reads
    bases = 0
    reads = 0

    # Iterate through each FASTQ file provided in the list
    for fp in files:
        # Open the compressed gzip file in text read mode ('rt')
        with gzip.open(fp, 'rt') as handle:
            # Loop through lines, starting the index at 1 for clean 4-line block matching
            for idx, line in enumerate(handle, 1):
                # FASTQ records have a strict 4-line structure: Header, Sequence, Spacer (+), Quality.
                # 'idx % 4 == 2' targets line 2 of every record, which is the actual DNA sequence.
                if idx % 4 == 2:
                    # Remove trailing newline characters and add the sequence length to the total
                    bases += len(line.rstrip())
                    # Increment the total read count
                    reads += 1

    # Initialize the target length with the pipeline's default fallback value
    trim_len = ${params.trimlength}
    
    # If reads were successfully found and parsed, calculate the true mathematical mean
    if reads > 0:
        trim_len = int(round(bases / reads))

    # Write the resulting calculated length to a local text file for provenance tracking
    with open('historical_trimlength.txt', 'w') as out:
        out.write(f"{trim_len}\\n")

    # Print to standard output so Nextflow can capture this value directly into a channel (.out)
    print(trim_len)
    PY
}

process BWA_MERGED {
    // map the merged reads
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq), val(mapper)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]
    
    if (mapper == "aln") {
        """
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${merged_fq} > ${sample_id}.sai
        
        bwa samse -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            ${params.reference} ${sample_id}.sai ${merged_fq} \
            | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    } else {
        """
        bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            -t ${task.cpus} ${params.reference} ${merged_fq} \
            | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    }
}

process BWA_UNMERGED {
    // map the unmerged reads
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2), val(mapper)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]

    if (mapper == "aln") {
        """
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r1} > ${sample_id}_R1.sai
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r2} > ${sample_id}_R2.sai

        bwa sampe \
            -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            ${params.reference} ${sample_id}_R1.sai ${sample_id}_R2.sai ${r1} ${r2} \
            | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    } else {
        """
        bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            -t ${task.cpus} ${params.reference} ${r1} ${r2} \
            | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    }
}

process MARKDUP_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_bam)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools view -@ ${task.cpus} -h ${merged_bam} | python3 ${params.rmdup_script} | samtools view -b -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.rmdup.bam
    """
}

process MARKDUP_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(unmerged_bam)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools collate -o ${sample_id}.sorted.namecollate.bam ${unmerged_bam}
    samtools fixmate -m ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam
    samtools sort -o ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.sorted.fixmate.bam
    samtools markdup -r ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.L${params.trimlength}.sorted.rmdup.bam
    
    # Cleanup
    rm ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam ${sample_id}.sorted.fixmate.positionsort.bam
    """
}

process MERGE_BAMS {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${params.trimlength}.bam"), path("${sample_name}.merged.L${params.trimlength}.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_name}.merged.L${params.trimlength}.bam ${bams}
    samtools index ${sample_name}.merged.L${params.trimlength}.bam
    """
}

process INDEL_REALN {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${params.trimlength}.realn.bam")

    script:
    """
    # Target Creator
    java -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator -R ${ref} -I ${bam} -o ${sample_name}.realn_targets.list -nt ${task.cpus}

    # Indel Realigner
    java -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner -R ${ref} -I ${bam} -targetIntervals ${sample_name}.realn_targets.list -o ${sample_name}.merged.L${params.trimlength}.realn.bam
    """
}

process INDEX_REALIGNED {
    tag "$sample_name"
    publishDir "${params.outdir}/data/bam", mode: 'copy'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}

process BAM_QC {
    tag "$sample_name"
    publishDir "${params.outdir}/results/stats", mode: 'copy'

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    path "${sample_name}.Q25.L${params.trimlength}.bam.dpstats.txt"

    script:
    """
    # Depth
    samtools depth -a -Q 30 -q 25 -b ${params.bed_file} ${bam} > ${sample_name}.Q25.bam.dp
    awk '{sum+=\$3} END { print sum/NR }' ${sample_name}.Q25.bam.dp > ${sample_name}.Q25.L${params.trimlength}.bam.dpstats.txt
    rm ${sample_name}.Q25.bam.dp
    """
}

process MAPDAMAGE {
    tag "$sample_name"
    publishDir "${params.outdir}/data/mapdamage", mode: 'copy', pattern: "*.{pdf,txt}"
    module 'container_env:mapdamage2' // not sure this is needed, but it doesn't hurt to load the module

    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref

    output:
    tuple val(sample_name), path("${sample_name}.rescaled.bam"), path("${sample_name}.rescaled.bam.bai"), emit: rescaled_indexed
    path "*.pdf", emit: plots, optional: true
    path "*.txt", emit: stats, optional: true

    script:
    """
    # Run mapDamage for damage assessment and rescaling
    crun mapDamage -i ${bam} -r ${ref} -d mapd_output -y
    
    # Copy rescaled BAM to standard output name
    if [ -f mapd_output/rescaled.bam ]; then
        cp mapd_output/rescaled.bam ${sample_name}.rescaled.bam
        samtools index ${sample_name}.rescaled.bam
    else
        echo "Error: mapDamage did not produce rescaled.bam" >&2
        exit 1
    fi
    
    # Copy report files
    cp mapd_output/*.pdf . 2>/dev/null || true
    cp mapd_output/*.txt . 2>/dev/null || true
    """
}

process AMBER_PREP {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("amber_input.txt"), path("${sample_name}.tags.MQ25.bam")

    script:
    """
    # Adding MD tags and filtering for MQ25
    samtools calmd -b ${bam} ${params.reference} > ${sample_name}.tags.bam
    samtools view -bq 25 ${sample_name}.tags.bam > ${sample_name}.tags.MQ25.bam
    
    # Create input file for AMBER
    echo -e "${sample_name}\t${sample_name}.tags.MQ25.bam" > amber_input.txt
    """
}

process AMBER {
    tag "$sample_name"
    publishDir "${params.outdir}/results/stats", mode: 'copy'

    input:
    tuple val(sample_name), path(amber_input), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}.amber_MQ25.pdf"), path("${sample_name}.amber_MQ25.txt")

    script:
    """
    # Run AMBER using the input file created in the previous step
    python3 ${params.amber_script} --bamfiles ${amber_input} --output ${sample_name}.amber_MQ25
    """
}

workflow MODERN_PIPELINE {
    take:
    modern_samples
    trim_length_ch
    mapper_choice_ch
    ref_ch
    ref_fai_ch
    ref_dict_ch

    main:
    modern_fastp = FASTP_MERGED(modern_samples)

    modern_merged_trim_input = modern_fastp.merged
        .combine(trim_length_ch)
        .map { sample_id, merged_fq, trim_len -> tuple(sample_id, merged_fq, trim_len) }
    modern_merged_split = SPLIT_MERGED(modern_merged_trim_input)

    modern_unmerged_fastp = FASTP_UNMERGED(modern_fastp.unmerged)
    modern_unmerged_trim_input = modern_unmerged_fastp[0]
        .combine(trim_length_ch)
        .map { sample_id, r1, r2, trim_len -> tuple(sample_id, r1, r2, trim_len) }
    modern_unmerged_trimmed = SEQTK_TRIM(modern_unmerged_trim_input)

    modern_merged_for_mapping = modern_merged_split
        .combine(mapper_choice_ch)
        .map { sample_id, merged_fq, mapper -> tuple(sample_id, merged_fq, mapper) }
    modern_unmerged_for_mapping = modern_unmerged_trimmed
        .combine(mapper_choice_ch)
        .map { sample_id, r1, r2, mapper -> tuple(sample_id, r1, r2, mapper) }

    QC_MERGED(modern_merged_split)
    QC_UNMERGED(modern_unmerged_trimmed)

    modern_bwa_merged = BWA_MERGED(modern_merged_for_mapping)
    modern_bwa_unmerged = BWA_UNMERGED(modern_unmerged_for_mapping)

    modern_markdup_merged = MARKDUP_MERGED(modern_bwa_merged)
    modern_markdup_unmerged = MARKDUP_UNMERGED(modern_bwa_unmerged)

    modern_markdup_merged
        .mix(modern_markdup_unmerged)
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            [sample_name, bam]
        }
        .groupTuple()
        .set { modern_bams_to_merge }
    modern_merge_bams = MERGE_BAMS(modern_bams_to_merge)

    modern_realn = INDEL_REALN(modern_merge_bams, ref_ch, ref_fai_ch, ref_dict_ch)
    modern_indexed = INDEX_REALIGNED(modern_realn)
    BAM_QC(modern_indexed)

    // Optional AMBER analysis
    if (params.run_modern_amber) {
        AMBER_PREP(modern_indexed)
        AMBER(AMBER_PREP.out)
    }

    emit:
    indexed = modern_indexed
}

workflow HISTORICAL_PIPELINE {
    // takes the output of FASTP directly, instead of the file paths like the modern pipeline
    // because the historical pipeline is designed to work with the cleaned reads from FASTP, 
    // which are already in the correct format for mapping and downstream processing.
    take:
    historical_merged_ch
    historical_unmerged_ch
    mapper_choice_ch
    ref_ch
    ref_fai_ch
    ref_dict_ch

    main:
    // redundant? since FASTP_MERGED already cleaned the reads with the same parameters
    historical_unmerged_fastp = FASTP_UNMERGED(historical_unmerged_ch)

    historical_merged_for_mapping = historical_merged_ch
        .combine(mapper_choice_ch)
        .map { sample_id, merged_fq, mapper -> tuple(sample_id, merged_fq, mapper) }
    historical_unmerged_for_mapping = historical_unmerged_fastp[0]
        .combine(mapper_choice_ch)
        .map { sample_id, r1, r2, mapper -> tuple(sample_id, r1, r2, mapper) }

    if (params.run_historical_fastqc) {
        QC_MERGED(historical_merged_ch)
        QC_UNMERGED(historical_unmerged_fastp[0])
    }

    historical_bwa_merged = BWA_MERGED(historical_merged_for_mapping)
    historical_bwa_unmerged = BWA_UNMERGED(historical_unmerged_for_mapping)

    historical_markdup_merged = MARKDUP_MERGED(historical_bwa_merged)
    historical_markdup_unmerged = MARKDUP_UNMERGED(historical_bwa_unmerged)

    historical_markdup_merged
        .mix(historical_markdup_unmerged)
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            [sample_name, bam]
        }
        .groupTuple()
        .set { historical_bams_to_merge }
    historical_merge_bams = MERGE_BAMS(historical_bams_to_merge)

    historical_realn = INDEL_REALN(historical_merge_bams, ref_ch, ref_fai_ch, ref_dict_ch)
    historical_indexed = INDEX_REALIGNED(historical_realn)
    BAM_QC(historical_indexed)

    historical_rescaled = Channel.empty()
    if (params.run_historical_mapdamage) {
        historical_mapdamage = MAPDAMAGE(historical_indexed, ref_ch)
        historical_rescaled = historical_mapdamage.rescaled_indexed
    }

    if (params.run_historical_amber) {
        AMBER_PREP(historical_indexed)
        AMBER(AMBER_PREP.out)
    }

    emit:
    indexed = historical_indexed
    rescaled = historical_rescaled
}

workflow {
    ref_ch      = Channel.fromPath(params.reference).first()
    ref_fai_ch  = Channel.fromPath("${params.reference}.fai").first()
    ref_dict_ch = Channel.fromPath(params.reference.replaceAll(/\.fasta$/, '.dict')).first()
    bed_ch      = Channel.fromPath(params.bed_file).first()

    // 0. Setup: Prepare reference repeat mask and calculate trim length from historical reads
    PREP_REFERENCE_REPEAT(ref_ch, bed_ch)

    // 1. Run FASTP on historical reads immediately to generate unmerged reads
    historical_fastp_input = historical_samples_ch
        .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }    
    historical_fastp = FASTP_MERGED(historical_fastp_input)

    // 2. Extract the paths from the unmerged output and pass to read length calculation
    // `.flatMap` breaks the tuples down into a flat stream, `.collect` packs them into a single list
    historical_unmerged_paths_ch = historical_fastp.unmerged
        .flatMap { sample_id, r1, r2 -> [r1, r2] }
        .collect()

    hist_trim_output = CALC_HISTORICAL_TRIMLEN(historical_unmerged_paths_ch)
    
    historical_trim_length_ch = hist_trim_output.trim_len
        .map { it.trim() }
        .map { it as Integer }

    // 3. Determine mapper based on trim length: aln for <= 80bp, mem for > 80bp
    // and log the output
    mapper_ch = historical_trim_length_ch.map { trim_len -> 
        def mapper = trim_len <= 80 ? "aln" : "mem"
        
        // This writes directly to the console output and the hidden .nextflow.log file
        log.info """
        ===================================================================
         PIPELINE CONFIGURATION RESOLVED:
         -> Calculated Historical Trim Length : ${trim_len} bp
         -> Chosen Alignment Mapper           : bwa ${mapper}
        ===================================================================
        """.stripIndent()

        return mapper
    }

    // 4. Fire off the pipelines
    modern_pipeline = MODERN_PIPELINE(modern_samples_ch, historical_trim_length_ch, mapper_ch, ref_ch, ref_fai_ch, ref_dict_ch)
    
    // Note: We pass the pre-calculated FASTP output directly into the Historical Pipeline
    historical_pipeline = HISTORICAL_PIPELINE(historical_fastp.merged, historical_fastp.unmerged, mapper_ch, ref_ch, ref_fai_ch, ref_dict_ch)

    // 5. Wire cleaned BAMs into ANGSD downstream modules
    modern_bams_for_angsd = modern_pipeline.indexed
        .map { sample_name, bam, bai -> tuple(sample_name, "modern", bam, bai) }

    if (params.run_historical_mapdamage && params.use_historical_rescaled) {
        historical_bams_for_angsd = historical_pipeline.rescaled
            .map { sample_name, bam, bai -> tuple(sample_name, "historical_rescaled", bam, bai) }
    } else {
        historical_bams_for_angsd = historical_pipeline.indexed
            .map { sample_name, bam, bai -> tuple(sample_name, "historical_original", bam, bai) }
    }

    bams_for_angsd = modern_bams_for_angsd.mix(historical_bams_for_angsd)
}