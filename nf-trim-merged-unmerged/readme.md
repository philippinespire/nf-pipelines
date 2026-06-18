## Overview

This pipeline trims fastq files and maps using ancient DNA-specific settings. It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster. The pipeline uses the same symbolic filename structure as the GenErode pipeline.

**Important for WAHAB HPC** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Installation

The pipeline can be installed directly from git

```bash

git clone https://github.com/mariannedehasque/nf-pipelines.git

```
This GitHub repository also contains the `nf-mapping-ancient-merged` and `nf-angsd-diversity` nextflow pipelines.

## Folder structure

## Parameters

### Inputfiles

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* Reference index files (.bwt, .ann, .sa, .pac, .ann, .amb )
* BED file with masked repeat regions
* File with the name of all fastq files to be processed 
* Directory containing all fastq files

All inputfiles can be copied or generated from the GenErode directory. Below is the code I used to generate the inputfiles.

```bash
# Activate bash shell

bash

# Change directory to the nextflow pipeline
cd ./nf-trim-merged-unmerged

# Create new directories
mkdir data
mkdir ./data/reference
mkdir ./data/symlinks
mkdir inputfiles

# Create softlinks to the raw fastq files
ln -s /Generode/data/raw_reads_symlinks/modern/*fastq.gz ./data/symlinks

# Create fastq filenames file. Manually adjust the file if necessary (e.g. if not all samples from GenErode are to be used)

ls ./data/symlinks/*fastq.gz | xargs -n1 basename | cut -d "_" -f1,2,3 | uniq > ./inputfiles/fastq_filenames.txt

# Create softlinks to reference, dict, and repma bed file
# Adjust the path to the reference if necessary

ln -s /Generode/reference/<reference>.fasta ./data/reference/
ln -s /Generode/reference/<reference>.dict ./data/reference/
ln -s /Generode/reference/<reference>.fasta.* ./data/reference/
ln -s /Generode/reference/<reference>.repma.bed ./data/reference/

```

**Important** The pipeline assumes that fastq files are named in the following way: 
`<sampleID>_<index>_<flowcellID>_R1.fastq.gz` eg. `TzoCMta031_1_22CVWFLT3_R1.fastq.gz`
Make sure that all fastq file names are unique and follow this structure. Otherwise the pipeline will fail.

### Configuration

In the `main.nf` file, adjust the parameters. If you used the same structure as outlined above, only the name of the reference should be added.

Params.trimlength will depend on the average read length of historical samples. I advice to keep the trim length < 89

```bash
// --- Default Parameters ---
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt"
params.indir        = "${projectDir}/data/symlinks"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/<reference>.fasta"
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed"
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/home/mdehasqu/TOOLS/AMBER/AMBER" // Ignore this.
params.bwa_threads  = 4
params.bam_q        = 1 // Mapping quality. Currently set to 1 simply to remove unmapped reads. 
params.trimlength   = 85
```

Since the `AMBER` tool requires special software and a unique conda environment on WAHAB, we will just skip this step for now. 
See below for instructions on how to run `AMBER` (not tested).

## To Run

I like to run nextflow pipelines in tmux so that I can keep it running in the background even when logging of.

To create a new tmux screen:

```bash
tmux new -s nextflow
```

!Good to know! Leaving a tmux sessions is notoriously difficult. To leave and reopen, type the following commands:

```bash
# To leave, literraly press these keys
Ctrl + B, followed by D

# To re-enter a session
tmux a -t nextflow

```

Then start the nextflow pipeline as follows (while in tmux):

```bash
# Load the nextflow module
module load container_env
module load nextflow

# Start the run
nextflow run main.nf -profile standard -resume

```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

## Including AMBER in the pipeline

1. Remove `//` in the AMBER step at the bottom of the `main.nf` file. The code should look like this.

```bash
    // Run AMBER (Prep -> Run)
    AMBER_PREP(INDEX_REALIGNED.out)
    AMBER(AMBER_PREP.out)

```

2. Update the `environment.yml` file

```bash
name: trimming
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - fastp=0.23.4
  - seqtk=1.4
  - bwa=0.7.17
  - samtools=1.19
  - fastqc=0.11.9
  - openjdk=11
  - fontconfig
  - freetype
  - python=3.9
  - numpy
  - matplotlib
  - pysam

```

3. Update the `nextflow.config` file and remove the reference to the conda environment.

```bash
// AMBER uses different conda environment
            withName: 'AMBER' {
                cpus = 2
                memory = '8 GB'
```