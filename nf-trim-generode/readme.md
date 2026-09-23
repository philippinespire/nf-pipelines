# nf-trim-generode: Historic & Modern DNA Mapping Pipeline

## Overview
This is a Nextflow DSL2 workflow designed for processing high-throughput sequencing data across two eras (historical and modern cohorts). It is an extension of `nf-trim-merged-unmerged` that incorporates GenErode features and takes both historical and modern reads. 

The pipeline performs:

* De novo repeat modeling, repeat masking, and CpG site identification on the reference genome.
* Generation of clean site-index files for downstream ANGSD analyses.
* Adapter trimming and overlapping read pair merging (fastp).
* Initial mapping of historical samples to calculate average mapped read length.
* Automated mapper selection (bwa aln vs bwa mem) based on empirical historical read lengths.
* Read trimming/splitting of modern reads down to historical read lengths to eliminate temporal length bias.
* Conditional re-mapping for historical samples if the optimal mapper differs from the first mapping.
* Custom duplicate removal, GATK Indel Realignment, depth statistics calculation, and optional FastQC, AMBER alignment evaluation, and mapDamage base quality rescaling.

The pipeline is configured to run on Old Dominion University's WAHAB cluster using Nextflow DSL2 and Slurm/Conda execution profiles.

---

## Installation & Environment Setup

The pipeline can be installed directly via git:

```bash
git clone https://github.com/philippinespire/nf-pipelines.git
cd ./nf-trim-merged-unmerged
```

> [!NOTE]
> This GitHub repository also contains other nextflow pipelines.

> [!IMPORTANT]
> Because you are likely to be cloning the nf-pipelines repo into your project repository, you should either add it to the `.gitignore` so it's not tracked, or delete the `nf-pipelines/.git` dir

### Personal Conda Prerequisite
> [!IMPORTANT]
> **For WAHAB HPC** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Folder structure
Before running, this is the setup:
```
nf-trim-merged-unmerged/  
├── data/  
│   ├── reference/         # Reference genome FASTA, .fai, .dict, and all BWA index layers  
│   └── symlinks/          # Symbolic links pointing to raw paired-end sequence files  
├── inputfiles/  
│   ├── samplesheet.csv    # Target metadata sheet (sample,era)  
│   └── fastq_filenames.txt # Fallback single-column sample ID list (optional)
├── scripts/  
│   ├── split_reads.sh     # Custom length splitting script for long, merged reads
│   └── samremovedup.py    # Custom python script for duplicate removal  
├── main.nf                # Master workflow execution script
├── mapping_modules.nf     # Module definitions for BWA ALN/MEM mapping
└── nextflow.config        # Process resource configs and cluster profiles
```

## Parameters & Configuration

All core parameters are declared in `main.nf` and can be customized in the script or overridden via the command line:

```
// --- Input & Path Parameters ---
params.samplesheet              = "${projectDir}/inputfiles/samplesheet.csv"           // Preferred way to provide sample metadata, including sample IDs and eras (modern or historical).
params.samples_file             = "${projectDir}/inputfiles/fastq_filenames.txt"       // Legacy way to provide sample metadata. One-column text file with sample IDs. All modern by default.
params.indir                    = "${projectDir}/data/symlinks"                        // Directory where raw FASTQ files are expected.
params.outdir                   = "${projectDir}/results"                              // Directory where all output files will be written.
params.reference                = "${projectDir}/data/reference/<reference>.fasta"     // Path to reference genome FASTA file.
params.bed_file                 = "${projectDir}/data/reference/<reference>.repma.bed" // Input bed file of repeats and CpG sites for downstream ANGSD analyses. Required if run_repeatmasking = false. Not used if run_repeatmasking = true.

// --- Pipeline Control & Logic Flags ---
params.historical_era           = "historical" // String identifier used in the samplesheet era column
params.historical_mapper        = "aln"        // Initial mapper choice for Pass 1 ("aln" or "mem")
params.run_repeatmasking        = true         // Set to false to bypass RepeatModeler/Masker and use params.bed_file
params.run_historical_fastqc    = true         // Run FastQC on historical reads after fastp
params.run_historical_mapdamage = true         // Run mapDamage assessment & base rescaling on historical BAMs
params.run_historical_amber     = true         // Run AMBER quality evaluation on historical BAMs
params.run_modern_amber         = true         // Run AMBER quality evaluation on modern BAMs

// --- Technical & Filtering Thresholds ---
params.bam_q                    = 25           // Mapping quality threshold for BAM filtering and length calculation
params.trimlength               = 85           // Fallback trim length if no historical reads pass mapping filters
params.bwa_threads              = 4            // Thread count per BWA process
```

Note that running repeat masking is a fairly slow process (a few hours). If you don't run it, however, you will need to specify a BED file.


### Inputfiles

The pipeline requires the following inputfiles:

* Reference (fasta format: .fasta, .fa, .fna)
* Reference index and dictionary files (.fai, .dict, .bwt, .pac, .ann, .amb, and .sa)
* CSV file with columns sample and era with the sample name and era (historical or modern) of all fastq files to be processed (`inputfiles/samplesheet.csv`). See [samplesheet_example.csv](examples/samplesheet_example.csv).
* Directory containing all fastq files

The pipeline assumes that fastq files are named in the following way: 
`<sampleID>_<index>_<flowcellID>_R1.fastq.gz` eg. `TzoCMta031_1_22CVWFLT3_R1.fastq.gz`  
and the paired `_R2.fastq.gz` file.  
Make sure that all fastq file names are unique and follow this structure. Otherwise the pipeline will fail.

Below is code that can help to generate the inputfiles.

```bash
# Activate bash shell. 
bash

# Change directory to the nextflow pipeline
cd ./nf-trim-generode

# Create new directories
mkdir data
mkdir ./data/reference
mkdir ./data/symlinks
mkdir inputfiles

# Create softlinks to the raw fastq files
ln -s /path/to/raw_reads/*fastq.gz ./data/symlinks

# Create samplesheet metadata file from the fastq files.
# Uses <sampleID>_<index>_<flowcellID> as the sample_id in the nextflow pipeline
# This assumes that the sample name in the fastq files have an A (Albatross) or C (contemporary) in the 4th position, e.g., TzoAMta031_1_22CVWFLT3 and TzoCMta031_1_22CVWFLT3
# Manually adjust the file if necessary (e.g. if not all samples are to be used)
(echo "sample,era"; ls ./data/symlinks/*fastq.gz | xargs -n1 basename | cut -d "_" -f1,2,3 | uniq | awk '{
    type = substr($0, 4, 1)
    if (type == "A") 
        print $0 ",historical"
    else if (type == "C") 
        print $0 ",modern"
    else 
        print $0 ",modern"
}') > ./inputfiles/samplesheet.csv

# Create a dictionary and index for the reference if it doesn't exist yet (adjust this to match your files)
salloc # the bwa in particular can be computationally intensive
module load container_env samtools/1.19
crun.samtools samtools dict <reference>.fna -o <reference>.dict
crun.samtools samtools faidx <reference>.fna
module unload container_env samtools/1.19
module load bwa
crun bwa index <reference>.fna
exit # leave the salloc if it was used


# Create softlinks to the genomic reference files:
# Adjust the path to the reference if necessary
ln -s /path/to/reference/<reference>.fasta ./data/reference/
ln -s /path/to/reference/<reference>.fasta.* ./data/reference/
```

## Pipeline Logic & Execution Flow
```
1. Reference Preparation (RepeatModeler -> RepeatMasker -> CpG Extraction -> ANGSD Sites Index)
                                  │
2. Historical Read Processing (fastp -> initial BWA Mapping -> Calculate Mapped Read Length)
                                  │
      ┌───────────────────────────┴───────────────────────────┐
      ▼                                                       ▼
Mapped Length ≤ 80 bp                                Mapped Length > 80 bp
  -> Select `bwa aln`                                  -> Select `bwa mem`
      │                                                       │
      └───────────────────────────┬───────────────────────────┘
                                  │
3. Modern Read Normalization (fastp -> Truncate/Split Reads to Historical Avg Length -> Mapping)
                                  │
4. BAM Post-Processing (Duplicate Removal -> GATK Indel Realignment -> Depth QC)
                                  │
5. Downstream Evaluation (Optional mapDamage Rescaling & AMBER Reports)
```

1. Reference Masking & ANGSD Site Indexing:
   * If `params.run_repeatmasking` = true, runs RepeatModeler for de novo repeat discovery, extracts CpG sites using a Python parser, and merges them via RepeatMasker.
   * Prepares 1-based inclusion sites files (`.repma.angsd.txt`), indexing them via angsd sites index for downstream genotype likelihood estimation.
2. Dynamic Mapper Selection & Two-Pass Mapping:
   * Historical reads undergo quality filtering and adapter trimming via fastp.
   * Reads map initially using `params.historical_mapper`. High-quality alignments ($Q \ge \text{params.bam\_q}$) are parsed to compute the true average mapped read length.
   * If average mapped read length is $\le 80\text{ bp}$, `bwa aln` is selected; if $> 80\text{ bp}$, `bwa mem` is selected.
   * If the chosen algorithm differs from the initial mapping, re-mapping is automatically executed.
3. Modern Read Normalization:
   * Modern long reads are truncated (`seqtk trimfq`) or split down the middle (`split_reads.sh`) to match the calculated historical average read length, eliminating temporal length bias signatures.
4. BAM Processing & Downstream QC:
   * Custom duplicate marking (`samremovedup.py` for merged reads, `samtools markdup` for unmerged pairs).
   * Local GATK Indel Realignment (`RealignerTargetCreator` and `IndelRealigner`).
   * Per-sample sequencing depth stats over target BED regions.
   * Optional `mapDamage` base-quality score rescaling for historical samples and `AMBER` alignment evaluation.


## To Run

It is helpful to run nextflow pipelines in tmux so that you can keep it running in the background even when logging off.

To create a new tmux screen:

```bash
tmux new -s nextflow
```

!Good to know! Leaving a tmux sessions is notoriously difficult. To leave and reopen, type the following commands:

```bash
# To leave, literally press these keys
Ctrl + B, followed by D

# To re-enter a session
tmux a -t nextflow

```

Then start the nextflow pipeline as follows (while in tmux):

```bash
# load bash with conda
bash
# Load the nextflow module
module load container_env nextflow

# Start the run
nextflow run main.nf -profile wahab -resume
```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

To override parameters directly from the command line:
```
nextflow run main.nf \
  --reference "./data/reference/my_genome.fasta" \
  --run_repeatmasking false \
  --bed_file "./data/reference/my_repeats.bed" \
  --historical_mapper "aln" \
  -profile wahab \
  -resume
```

## Output
Outputs are written to `params.outdir` (default: `./results`):
```
results/
├── amber/               # AMBER alignment quality plots (.pdf) and summary stats (.txt)
├── data/
│   ├── bam/             # Final sorted, deduplicated, and indel-realigned BAM and BAI files
│   ├── bam_rescaled/    # Post-mapDamage base-quality rescaled BAM and BAI files (historical)
│   ├── fastq/           # Length-trimmed/split modern FASTQ files
│   └── reference/       # Masking tracks (.combined_mask.bed) and ANGSD site index files:
│                        #   ├── <ref>.repma.angsd.txt (.idx, .bin)
│                        #   ├── <ref>.regions
│                        #   └── <ref>.chrs
│       └── modeler/     # De novo repeat library (consensi.fa, families.stk)
├── depth/               # Average coverage depth stats over target regions (.dpstats.txt)
├── fastp/               # Fastp HTML and JSON trimming reports
├── fastqc/              # FastQC reports (.html, .zip) for trimmed reads
├── mapdamage/           # mapDamage postmortem damage plots (.pdf) and statistic logs (.txt)
└── stats/               # Calculated historical average read length (historical_trimlength.txt)
```

## Software Stack
The pipeline uses:

* fastp (v0.23.4) — Quality-trimming and overlapping mate pairing/merging.
* seqtk (v1.4) — Down-sampling and structural end-read truncation mapping.
* bwa (v0.7.17) — Dual aln/mem short-read genome matching engine.
* samtools (v1.19) — Indexing, sorting, track filtering, and deep map coverage summaries.
* GATK3 — Target-based local indel realignment (RealignerTargetCreator / IndelRealigner).
* RepeatModeler / RepeatMasker — Genomic repeat identification and masking.
* fastqc (v0.11.9) — Per-pass data quality assurance metrics.
* mapDamage (v2.2) — Postmortem historical damage assessment and base-quality score rescaling.
* AMBER (v2.0) — Target alignment quality extraction evaluations.
* angsd (v0.94) — Multi-individual genotype processing database assembly.
