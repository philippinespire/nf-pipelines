// Marianne Dehasque and Malin Pinsky, 2026
// Assistance from GPT5.3 and and Gemini.

nextflow.enable.dsl=2
include { BWA_MERGED as BWA_MERGED_PASS1; BWA_MERGED as BWA_MERGED_PASS2 } from './mapping_modules.nf'
include { BWA_UNMERGED as BWA_UNMERGED_PASS1; BWA_UNMERGED as BWA_UNMERGED_PASS2 } from './mapping_modules.nf'

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // This is the preferred way to provide sample metadata, including sample IDs and eras (modern or historical). If this file is not found, the pipeline will fall back to using the legacy samples_file parameter.
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt" // This is the legacy way to provide sample metadata. It should be a one-column text file with sample IDs. All samples will be treated as modern if this file is used.
params.indir        = "${projectDir}/data/symlinks" // This is the directory where the raw FASTQ files are expected to be located. The pipeline will look for files named <sample_id>_R1.fastq.gz and <sample_id>_R2.fastq.gz in this directory.
params.outdir       = "${projectDir}/results" // This is the directory where all output files will be written. The pipeline will create subdirectories for different types of output (e.g., fastq, bam, stats).
params.reference    = "${projectDir}/data/reference/<reference>.fasta" // This is the path to the reference genome FASTA file that will be used for read mapping. The <reference> placeholder should be replaced with the actual reference name (e.g., hg19, mm10).
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed" // This is the path to the BED file containing repeat-masked regions of the reference genome. The <reference> placeholder should be replaced with the actual reference name (e.g., hg19, mm10). This file is used for filtering reads during mapping and for calculating depth statistics.
params.reference_prefix         = params.reference.tokenize('/').last().replaceAll(/\.(fa|fasta|fna)$/, '') // This extracts the base name of the reference file without the directory path and without the file extension. It is used for naming output files related to the reference genome.
params.historical_era           = "historical"
params.historical_mapper        = "aln" // Default starting point for historical read mapping
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

process CALC_HISTORICAL_TRIMLEN_BAM {
    publishDir "${params.outdir}/data/stats", mode: 'copy'
    tag "Calculating true mapped length"

    input:
    path bams // Nextflow stages all collected BAMs into the working directory

    output:
    stdout emit: trim_len
    path "historical_trimlength.txt", emit: trim_file

    script:
    """
    # Loop through all BAMs and stream the mapped reads into a single awk process
    for bam in ${bams}; do
        samtools view -F 2308 "\$bam"
    done | awk '
      BEGIN { bases=0; reads=0 }
      {
        # Column 10 is the sequence. Ignore if it is a missing "*"
        if (\$10 != "*") {
          bases += length(\$10)
          reads++
        }
      }
      END {
        if (reads > 0) {
          # Calculate average and round to nearest integer
          printf "%d\\n", bases/reads + 0.5
        } else {
          # Fallback to default if no mapped reads exist
          print ${params.trimlength}
        }
      }' | tee historical_trimlength.txt
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

process MAPDAMAGE {
    tag "$sample_name"
    publishDir "${params.outdir}/data/mapdamage", mode: 'copy'
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
    crun mapDamage -i ${bam} -r ${ref} -d mapd_output --rescale --merge-reference-sequences
    
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
    // this takes raw reads as input
    take:
    modern_samples
    trim_length_ch
    mapper_choice_ch
    ref_ch
    ref_fai_ch
    ref_dict_ch

    main:
    // 1. Run the initial FASTP_MERGED
    modern_fastp = FASTP_MERGED(modern_samples)

    // 2. Handle Merged Reads (Split them if they are longer than the trim length)
    modern_merged_trim_input = modern_fastp.merged
        .combine(trim_length_ch)
        .map { sample_id, merged_fq, trim_len -> tuple(sample_id, merged_fq, trim_len) }
    modern_merged_split = SPLIT_MERGED(modern_merged_trim_input)

    // 3. Handle Unmerged Reads (Truncate them directly with SEQTK)
    modern_unmerged_trim_input = modern_fastp.unmerged
        .combine(trim_length_ch)
        .map { sample_id, r1, r2, trim_len -> tuple(sample_id, r1, r2, trim_len) }
    
    modern_unmerged_trimmed = SEQTK_TRIM(modern_unmerged_trim_input)

    // 4. Combine with the mapper choice
    modern_merged_for_mapping = modern_merged_split
        .combine(mapper_choice_ch)
        .map { sample_id, merged_fq, mapper -> tuple(sample_id, merged_fq, mapper) }
    modern_unmerged_for_mapping = modern_unmerged_trimmed
        .combine(mapper_choice_ch)
        .map { sample_id, r1, r2, mapper -> tuple(sample_id, r1, r2, mapper) }

    // 5. QC steps
    QC_MERGED(modern_merged_split)
    QC_UNMERGED(modern_unmerged_trimmed)

    // 6. Mapping
    modern_bwa_merged = BWA_MERGED_PASS1(modern_merged_for_mapping)
    modern_bwa_unmerged = BWA_UNMERGED_PASS1(modern_unmerged_for_mapping)

    // 7. Mark duplicates
    modern_markdup_merged = MARKDUP_MERGED(modern_bwa_merged)
    modern_markdup_unmerged = MARKDUP_UNMERGED(modern_bwa_unmerged)

    // 8. Merge BAMs
    modern_markdup_merged
        .mix(modern_markdup_unmerged)
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            [sample_name, bam]
        }
        .groupTuple()
        .set { modern_bams_to_merge }
    modern_merge_bams = MERGE_BAMS(modern_bams_to_merge)

    // 9. Indel realignment and indexing
    modern_realn = INDEL_REALN(modern_merge_bams, ref_ch, ref_fai_ch, ref_dict_ch)
    modern_indexed = INDEX_REALIGNED(modern_realn)
    BAM_QC(modern_indexed)

    // 10. Optional AMBER analysis
    if (params.run_modern_amber) {
        AMBER_PREP(modern_indexed)
        AMBER(AMBER_PREP.out)
    }

    emit:
    indexed = modern_indexed
}

workflow HISTORICAL_PIPELINE {
    // This takes mapped reads as input. The QC and mapping is done in the main workflow.
    take:
    historical_bwa_merged
    historical_bwa_unmerged
    ref_ch
    ref_fai_ch
    ref_dict_ch

    main:
    // 1. Mark duplicates
    historical_markdup_merged = MARKDUP_MERGED(historical_bwa_merged)
    historical_markdup_unmerged = MARKDUP_UNMERGED(historical_bwa_unmerged)

    // 2. Merge BAMs
    historical_markdup_merged
        .mix(historical_markdup_unmerged)
        .map { id, bam -> [id.split('_')[0], bam] }
        .groupTuple()
        .set { historical_bams_to_merge }
        
    historical_merge_bams = MERGE_BAMS(historical_bams_to_merge)

    // 3. Indel realignment and indexing
    historical_realn = INDEL_REALN(historical_merge_bams, ref_ch, ref_fai_ch, ref_dict_ch)
    historical_indexed = INDEX_REALIGNED(historical_realn)
    BAM_QC(historical_indexed)

    // 4. Optional MapDamage
    historical_rescaled = Channel.empty()
    if (params.run_historical_mapdamage) {
        historical_mapdamage = MAPDAMAGE(historical_indexed, ref_ch)
        historical_rescaled = historical_mapdamage.rescaled_indexed
    }

    // 5. Option AMBER analysis
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

    // 1. FASTP on Historical Reads
    historical_fastp_input = historical_samples_ch.map { id, r1, r2 -> tuple(id, [r1, r2]) }    
    historical_fastp = FASTP_MERGED(historical_fastp_input)
    
    // 2. Option QC on historical reads
    if (params.run_historical_fastqc) {
        QC_MERGED(historical_fastp.merged)
        QC_UNMERGED(historical_fastp.unmerged)
    }

    // 3. Initial historical mapping based on user-chosen mapper
    pass1_mapper_ch = Channel.value(params.historical_mapper)
    
    hist_merged_pass1_in = historical_fastp.merged.combine(pass1_mapper_ch)
    hist_unmerged_pass1_in = historical_fastp.unmerged.combine(pass1_mapper_ch)

    hist_bwa_merged_pass1 = BWA_MERGED_PASS1(hist_merged_pass1_in)
    hist_bwa_unmerged_pass1 = BWA_UNMERGED_PASS1(hist_unmerged_pass1_in)

    // 3.1 Gather BAMs for calculation
    pass1_bams_ch = hist_bwa_merged_pass1.map { id, bam -> bam }
        .mix(hist_bwa_unmerged_pass1.map { id, bam -> bam })
        .collect()

    // 4. Calculate mapped length and choose mapper
    hist_trim_output = CALC_HISTORICAL_TRIMLEN_BAM(pass1_bams_ch)
    historical_trim_length_ch = hist_trim_output.trim_len.map { it.trim() as Integer }

    mapper_ch = historical_trim_length_ch.map { trim_len -> 
        def chosen = trim_len <= 80 ? "aln" : "mem"
        log.info """
        ===================================================================
         PASS 1 MAPPING COMPLETE:
         -> User-Specified Mapper             : bwa ${params.historical_mapper}
         -> Calculated Mapped Read Length     : ${trim_len} bp
         -> Optimal Chosen Mapper             : bwa ${chosen}
        ===================================================================
        """.stripIndent()
        return chosen
    }

    // 5. Condition remapping of historical reads    
    // 5.1. Gates: These filters only let data pass through if the mappers don't match
    hist_merged_pass2_in = historical_fastp.merged
        .combine(mapper_ch)
        .filter { id, fq, chosen -> chosen != params.historical_mapper }
        
    hist_unmerged_pass2_in = historical_fastp.unmerged
        .combine(mapper_ch)
        .filter { id, r1, r2, chosen -> chosen != params.historical_mapper }

    // 5.2 Run mapping a second time (Nextflow will just ignore this if the channels above are empty!)
    hist_bwa_merged_pass2 = BWA_MERGED_PASS2(hist_merged_pass2_in)
    hist_bwa_unmerged_pass2 = BWA_UNMERGED_PASS2(hist_unmerged_pass2_in)

    // 5.3 Merge logic: Route the correct BAMs to the final downstream pipeline
    final_hist_bwa_merged = hist_bwa_merged_pass1
        .combine(mapper_ch)
        .filter { id, bam, chosen -> chosen == params.historical_mapper } // Keep Pass 1 if mappers match
        .map { id, bam, chosen -> tuple(id, bam) }
        .mix(hist_bwa_merged_pass2) // Mix in Pass 2 (if it ran)

    final_hist_bwa_unmerged = hist_bwa_unmerged_pass1
        .combine(mapper_ch)
        .filter { id, bam, chosen -> chosen == params.historical_mapper } // Keep Pass 1 if mappers match
        .map { id, bam, chosen -> tuple(id, bam) }
        .mix(hist_bwa_unmerged_pass2) // Mix in Pass 2 (if it ran)
    
    // 6. Modern pipeline uses the new calculated length and chosen mapper
    modern_pipeline = MODERN_PIPELINE(modern_samples_ch, historical_trim_length_ch, mapper_ch, ref_ch, ref_fai_ch, ref_dict_ch)
    
    // 7. Historical pipeline takes the finalized BAMs directly
    historical_pipeline = HISTORICAL_PIPELINE(final_hist_bwa_merged, final_hist_bwa_unmerged, ref_ch, ref_fai_ch, ref_dict_ch)
    
    // 8. Wire cleaned BAMs into ANGSD downstream modules
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