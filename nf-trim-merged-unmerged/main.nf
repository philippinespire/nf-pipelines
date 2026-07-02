nextflow.enable.dsl=2

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt"
params.indir        = "${projectDir}/data/symlinks"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/<reference>.fasta"
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed"
params.reference_prefix = params.reference.tokenize('/').last().replaceAll(/\.(fa|fasta|fna)$/, '')
params.historical_era = "historical"
params.run_historical_fastqc = false
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/home/mdehasqu/TOOLS/AMBER/AMBER" // Ignore this.
params.bwa_threads  = 4
params.bam_q        = 1 // Mapping quality. Currently set to 1 simply to remove unmapped reads. 
params.trimlength   = 85

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

sample_metadata_ch.into { sample_metadata_for_modern; sample_metadata_for_historical }

modern_samples_ch = sample_metadata_for_modern
    .filter { sample_id, era -> era != params.historical_era }
    .map { sample_id, era -> tuple(sample_id, resolve_reads(sample_id)) }

historical_samples_ch = sample_metadata_for_historical
    .filter { sample_id, era -> era == params.historical_era }
    .map { sample_id, era ->
        def reads = resolve_reads(sample_id)
        tuple(sample_id, reads[0], reads[1])
    }

historical_samples_ch.into { historical_samples_for_qc; historical_samples_for_mapping; historical_samples_for_trimlen }

historical_read_paths_ch = historical_samples_for_trimlen
    .flatMap { sample_id, r1, r2 -> [r1.toString(), r2.toString()] }
    .collect()

// --- Processes ---

process FASTP_MERGED {
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
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq), val(trim_len)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${trim_len}.fastq.gz")

    script:
    """
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

files = [${reads_list}]
bases = 0
reads = 0

for fp in files:
    with gzip.open(fp, 'rt') as handle:
        for idx, line in enumerate(handle, 1):
            if idx % 4 == 2:
                bases += len(line.rstrip())
                reads += 1

trim_len = ${params.trimlength}
if reads > 0:
    trim_len = int(round(bases / reads))

with open('historical_trimlength.txt', 'w') as out:
    out.write(f"{trim_len}\\n")

print(trim_len)
PY
    """
}

process BWA_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]
    
    """
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${merged_fq} > ${sample_id}.sai
    
    bwa samse -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${sample_id}.sai ${merged_fq} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process BWA_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]

    """
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r1} > ${sample_id}_R1.sai
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r2} > ${sample_id}_R2.sai

    bwa sampe \
        -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${sample_id}_R1.sai ${sample_id}_R2.sai ${r1} ${r2} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
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

workflow {
    ref_ch      = Channel.fromPath(params.reference).first()
    ref_fai_ch  = Channel.fromPath("${params.reference}.fai").first()
    ref_dict_ch = Channel.fromPath(params.reference.replaceAll(/\.fasta$/, '.dict')).first()
    bed_ch      = Channel.fromPath(params.bed_file).first()

    PREP_REFERENCE_REPEAT(ref_ch, bed_ch)
    historical_trim_length_ch = CALC_HISTORICAL_TRIMLEN(historical_read_paths_ch)
        .out
        .trim_len
        .map { it.trim() as Integer }

    // 1. Modern branch: trim/merge/split while keeping merged and unmerged streams separate.
    modern_fastp = FASTP_MERGED(modern_samples_ch)
    modern_merged_trim_input = modern_fastp.out.merged
        .combine(historical_trim_length_ch)
        .map { sample_id, merged_fq, trim_len -> tuple(sample_id, merged_fq, trim_len) }
    modern_merged_split = SPLIT_MERGED(modern_merged_trim_input)

    modern_unmerged_fastp = FASTP_UNMERGED(modern_fastp.out.unmerged)
    modern_unmerged_trim_input = modern_unmerged_fastp.out[0]
        .combine(historical_trim_length_ch)
        .map { sample_id, r1, r2, trim_len -> tuple(sample_id, r1, r2, trim_len) }
    modern_unmerged_trimmed = SEQTK_TRIM(modern_unmerged_trim_input)
    
    // 2. QC
    QC_MERGED(modern_merged_split.out)
    QC_UNMERGED(modern_unmerged_trimmed.out)
    if (params.run_historical_fastqc) {
        QC_UNMERGED(historical_samples_for_qc)
    }
    
    // 3. Mapping (Per lane), still separated by stream type.
    modern_bwa_merged = BWA_MERGED(modern_merged_split.out)
    modern_bwa_unmerged = BWA_UNMERGED(modern_unmerged_trimmed.out)
    historical_bwa_unmerged = BWA_UNMERGED(historical_samples_for_mapping)

    // 4. Mark Duplicates (Split Processes)
    modern_markdup_merged = MARKDUP_MERGED(modern_bwa_merged.out)
    modern_markdup_unmerged = MARKDUP_UNMERGED(modern_bwa_unmerged.out)
    historical_markdup_unmerged = MARKDUP_UNMERGED(historical_bwa_unmerged.out)

    // 5. Merge BAMs (Per Biological Sample)
    // Mix both streams (merged & unmerged), convert Lane ID to Sample Name, then Group
    modern_markdup_merged.out
        .mix(modern_markdup_unmerged.out)
        .map { id, bam -> 
            def sample_name = id.split('_')[0] // e.g. "TzoCMta031" from "TzoCMta031_1_22CVWFLT3L3"
            return [ sample_name, bam ]
        }
        .groupTuple() // Groups all BAMs (merged & unmerged) for "TzoCMta031"
        .set { bams_to_merge }

    modern_merge_bams = MERGE_BAMS(bams_to_merge)

    historical_bams_to_merge = historical_markdup_unmerged.out
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            return [ sample_name, bam ]
        }
        .groupTuple()

    historical_merge_bams = MERGE_BAMS(historical_bams_to_merge)
    
    // 6. Indel Realign & QC
    modern_realn = INDEL_REALN(modern_merge_bams.out, ref_ch, ref_fai_ch, ref_dict_ch)
    modern_indexed = INDEX_REALIGNED(modern_realn.out)
    
    // Run Depth QC
    BAM_QC(modern_indexed.out)

    historical_realn = INDEL_REALN(historical_merge_bams.out, ref_ch, ref_fai_ch, ref_dict_ch)
    historical_indexed = INDEX_REALIGNED(historical_realn.out)
    BAM_QC(historical_indexed.out)

    // Run AMBER (Prep -> Run)
    // AMBER_PREP(INDEX_REALIGNED.out)
    // AMBER(AMBER_PREP.out)
}