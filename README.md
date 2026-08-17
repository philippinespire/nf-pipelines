## Overview
Repository to document custom nextflow pipelines. Pipelines were written for the WAHAB cluster.

Documentation and testing is work in progress.

## Content

```
nf-pipelines/
├── README.md
├── docs/
│   ├── lab_notebook.md
│   └── nf-trim-generode-porting-plan.md # some notes on creating nf-trim-generode
├── environments/
│   └── nf-angsd-diversity.yml  # Conda environment file for Wahab used by some pipelines
├── nf-trim-merged-unmerged/    # Takes modern reads. Trims modern to an average length 
│                               # and maps to the reference. Usually run after Generode.
├── nf-trim-bwamem/             # Like nf-trim-merged-unmerged, but uses bwamem
├── nf-trim-generode/           # Takes historical and modern reads. Trims modern to average
│                               # historical length and maps both with the same mapper. Also
│                               # masks repeats in the reference. 
│                               # A merger of Generode and nf-trim-merged-unmerged. 
│                               # Runs QA/QC with fastp, fastqc, amber, and mapdamage 
├── nf-angsd-diversity/         # Runs angsd for genotypes, calculates diversity on snps and a
│                               # PCA. Run this after nf-trim-*
└── nf-angsd-selection/         # Like nf-angsd-diversity. Adds a selection scan, 
                                # ld-pruning, admixture and fst. Runs diversity on all 
                                # putatively neutral sites (including monomorphic)

```
## QC steps before running the pipelines
The original pipelines (nf-trim-merged-unmerged and nf-angsd-diversity) were developed as follow up steps for GenErode. Later pipelines (nf-trim-generode) was developed to avoid the need for running Generode first. See the `readme.md` files under each pipeline for more information on using each one.

Before running the original pipelines (e.g., `nf-trim-merged-unmerged` and `nf-angsd-diversity`), it is advised to run the following QC scripts on the GenErode output to estimate average readlength, check sample (breadth of) coverage and mapping bias. These are also used as input for the various pipelines. I also provided a script to generate publication-ready sequencing statistics. Note that all these scripts were originally written by and for Marianne on the Old Dominion WAHAB cluster.

The different scripts are:
- Average readlength
- Average genome-wide coverage
- [AMBER](https://github.com/tvandervalk/AMBER) plots
- Sequencing efficiency

### Average readlength

```bash

docker=/cm/shared/containers/docker
samtools_image=$docker/biocontainers-samtools:v1.9-4-deb_cv1.sif

DIR=/path/to/output/directory
INPUT_BAM=/path/to/bamfile
NAME=<sample ID>
BED=/path/to/reference.repma.bed

echo $NAME

mkdir -p $DIR/results

date

echo "Calculating depth..."

singularity exec $samtools_image samtools depth -a -Q 30 -q 25 $INPUT_BAM > $DIR/results/${NAME}.Q25.bam.dp
awk '{{sum+=$3}} END {{ print sum/NR }}' $DIR/results/${NAME}.Q25.bam.dp > $DIR/results/${NAME}.Q25.bam.dpstats.txt

```
### Average genome-wide coverage

```bash

docker=/cm/shared/containers/docker
samtools_image=$docker/biocontainers-samtools:v1.9-4-deb_cv1.sif

DIR=/path/to/output/directory
INPUT_BAM=/path/to/bamfile
NAME=<sample ID>
BED=/path/to/reference.repma.bed

echo $NAME

mkdir -p $DIR/results

echo "Calculating average read length..."

singularity exec $samtools_image samtools view -h -q 25 -L $BED $INPUT_BAM	| grep -v '@' | awk '{print length($10)}'|  awk '{ total += $1; count++ } END { print total/count }' > $DIR/results/${NAME}.Q25.avRL.txt

```
### AMBER plots

#### Create conda environment (only first time)

```bash
#Activate bash
bash

#Create AMBER environment
conda create -n amber -c conda-forge -c bioconda pysam matplotlib numpy
```

>Important for WAHAB HPC: AMBER installation and activation is easiest with a personal conda environment. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

#### Running AMBER

```bash
#!/bin/bash -le
#SBATCH --job-name=AMBER
#SBATCH -o bamqc-%A_%a.out
#SBATCH --cpus-per-task=4

source ~/.bashrc
conda activate AMBER_environment

AMBER=/archive/carpenterlab/pire/softwares/AMBER/AMBER
REFERENCE=/path/to/reference.fasta

echo "Running AMBER..."

date

#Adding MD tags and filtering for MQ25
singularity exec $samtools_image samtools calmd -b $INPUT_BAM $REFERENCE > $DIR/data/$NAME.merged.rmdup.merged.realn.tags.bam
singularity exec $samtools_image samtools view -bq 25 $DIR/data/$NAME.merged.rmdup.merged.realn.tags.bam > $DIR/data/$NAME.merged.rmdup.merged.realn.tags.MQ25.bam

#Create input file for AMBER
echo -e "${NAME}\t$DIR/data/${NAME}.merged.rmdup.merged.realn.tags.MQ25.bam" > $DIR/data/${NAME}.amber_input_MQ25.txt

#run AMBER
crun.conda -p ~/envs/AMBER python3 $AMBER --bamfiles $DIR/data/${NAME}.amber_input_MQ25.txt --output $DIR/results/$NAME.amber_MQ25

# Cleanup
rm $DIR/data/${NAME}.amber_input_MQ25.txt
rm $DIR/data/$NAME.merged.rmdup.merged.realn.tags.MQ25.bam
rm $DIR/data/$NAME.merged.rmdup.merged.realn.tags.bam
```
### Sequencing efficiency

The scripts below are used to calculate endogenous content, complexity, coverage,and sequencing efficiency.
These metrics are important to report as metadata in publications and to determine sequencing effort after screening runs.

- collect_stats.sh
- collect_stats_modern.sh
- calculate.awk (called by collect_stats scripts)

**Example output**

| SAMPLENAME | total_seqs | total_mapped | total_uniq | total_MQ25 | endogenous | complexity | grr | total_cov | read_min | read_max | read_median | read_mean |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SdA0605708AWGSZ | 458713 | 69361 | 64173 | 50365 | 0.15121 | 0.92520 | 0.10980 | 0.02332 | 31 | 269 | 81 | 87.31 |
| SdA0609312EWGSZ | 5.446e+06 | 791050 | 675387 | 588302 | 0.14525 | 0.85379 | 0.10802 | 0.32385 | 31 | 269 | 96 | 102.33 |
| SdA0608411DWGSZ | 7.293e+06 | 361942 | 303596 | 259758 | 0.04963 | 0.83880 | 0.03562 | 0.15671 | 31 | 269 | 106 | 112.01 |

