## Overview of nf-angsd-selection

This pipeline runs (optional) chi-squared selection scans from data across two time-points, does (optional) ld-pruning, calculates diversity estimates within the ANGSD framework on all putatively neutral sites (including monomorphic sites), and then makes PCA and admixture plots and calculates Fst with ANGSD on putatively neutral ld-pruned sites. It only considers sites that meet minimum individual and maximum depth limits in each of the population-eras (set as a fraction and a multiplier, respectively).

It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster (use `-profile wahab`). The pipeline starts from BAM files, such as those output by nf-trim-generode. It was developed from nf-angsd-diversity.

**Important for WAHAB HPC:** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Installation

The pipeline can be installed along with the others in this repo directly from git

```bash

git clone https://github.com/mariannedehasque/nf-pipelines.git

```

## Required files

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* BED file with masked repeat regions (eg, output by Generode or by nf-trim-generode)
* BAM files

In addition, we need to create the following files too (see instructions below):
* Samplesheet
* ANGSD sites file
* Bam inputfile
* Contig file

### Samplesheet

This is a samplesheet unique to this pipeline. For an example, see ./inputfiles/samplesheet.csv

The samplesheet is a comma-separated file (csv) and consists of the following columns:

| sample     | bam                                               | pop     | era    | region    |
|------------|---------------------------------------------------|---------|--------|-----------|
| TzoCMal001 | /path/to/bam/TzoCMal001.merged.realn.bam          | TzoCMal | modern | Malampaya |
| TzoABol002 | /path/to/bam/TzoABol002.merged.realn.bam          | TzoABol | historic | Malampaya |
| TzoCMal003 | /path/to/bam/TzoCMal003.merged.realn.bam          | TzoCMal | modern | Malampaya |

Important:
* `era` takes one of two arguments: modern or historic
* `region` is used to define which historic and modern sites should be directly compared. Use the same value here if you want to directly compare two sites.
* `pop` is used to define the populations for which to calculate diversity matrixes

If you want to use mapdamage-rescaled BAM files for historical samples (eg, from Generode or nf-trim-generode), point to those files in this input file.


### ANGSD sites file
The ANGSD sites file contains information on the sites that will be analyzed. We will use this file to remove repeat regions from the analysis. This files corresponds to the `-sites` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Sites).

To create the sites file from the GenErode output (after loading all modules):

```bash
awk '{print $1"\t"$2+1"\t"$3}' ./path/to/reference/reference.repma.bed > ./path/to/reference/reference.repma.angsd.txt

angsd sites index ./path/to/reference/reference.repma.angsd.txt
```

Alternatively, this file is created automatically as part of the nf-trim-generode pipeline and output in `results/data/reference/*.repma.angsd.txt`.

### Bam inputfile
The filelist is a file containing the full path for each bam file with one filename per row.

This corresponds to the input given with the `-bam` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input).

If you want to use mapdamage-rescaled BAM files for historical samples (eg, from Generode or nf-trim-generode), point to those files in this input file.

### Contig file
Specify the contigs/regions for which to run the pipeline. One contig/region per line. This corresponds to the `-r` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input). The pipeline will submit a job per contig/region.

## Configuration

In the `main.nf` file, adjust the file paths. Make sure all files are in the correct directory.

In addition, the following parameters should be defined:
* min_ind_frac: The fraction of individuals that have to be represented at a locus/site to use it. This is applied per population-era. I use 70-80% number of total individuals here, though this also depends on the dataset.
* max_depth_mult: The multiple of the average depth to use as a ceiling in order to include a locus/site. This is applied per population-era. I use 10X the expected total coverage here.
* ld_prune: Set to true if you want to prune out loci in linkage disequilibrium (often a good idea). The PCA, admixture, and diversity calculations will then use the LD-pruned set of loci.
* max_kb_dist: LD pruning uses a sliding window to test for disequilibrium. This sets the window width. 50 kb is usually good. Narrower will run faster, wider will run more slowly.
* min_weight: Loci that are correlated (r2) more than this threshold will be trimmed out. 0.2 is often a good value.
* run_selection: Set to true to run chi-squared selection scan
* generations: Number of generations to use for Ne calculations. Default is 114 for Albatross to contemporary populations with a 1 year generation time. Adjust for other species as needed.
* fdr_cutoff: False detection rate (FDR) cutoff for iterative selection scan (default 0.05)
* max_rounds: Maximum number of rounds to run for iterative selection scan (default 20)
* n_boot: Number of bootstrap replicates to run for Ne calculations (default 1000)
* min_ind: Minimum number of individuals to include in Ne calculations (default 4)


```bash
// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/data/reference/reference.contigs.txt"
params.outdir      = "${projectDir}/results"
params.reference   = "${projectDir}/data/reference/reference.fasta"
params.bed_file    = "${projectDir}/data/reference/reference.repma.angsd.txt"
params.species     = "Sor"
params.max_depth_mult = 10  // The depth multipler
params.min_ind_frac   = 0.7 // The minimum fraction of individuals
params.ld_prune    = true   // Set to true to enable LD pruning
params.max_kb_dist = 50     // Maximum pairwise distance in kb to test for LD if pruning
params.min_weight  = 0.2    // Minimum r2 threshold for pruning filter
params.run_selection = true // Set to true to run ACER select scan
params.generations = 114    // Number of generations to use for Ne calculations. Default is 114 for Albatross to contemporary populations with a 1 year generation time. Adjust for other species as needed.
params.fdr_cutoff  = 0.05   // FDR cutoff for iterative selection scan (default 0.05)
params.max_rounds  = 20     // Maximum number of rounds to run for iterative selection scan (default 20)
params.n_boot      = 1000   // Number of bootstrap replicates to run for Ne calculations (default 1000)
params.min_ind     = 4      // Minimum number of individuals to include in Ne calculations

```

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

# Activate bash
bash

# Start the run
nextflow run main.nf -profile standard -resume

```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

## Output
## Pipeline Output Directory Structure

All pipeline outputs are written to the directory specified by `--outdir` (default: `./results/`). Large binary/intermediate files are isolated in the `large_data/` directory so they can be easily managed or excluded via `.gitignore`.

```text
results/
├── inputfiles/                # Input BAM file lists
├── sites/                     # Genomic site position sets
├── selection/                 # ACER iterative selection scan outputs
├── pcangsd/                   # Population structure (PCA & Admixture)
├── diversity/                 # Summary nucleotide diversity statistics
├── fst/                       # FST values and 2D SFS
└── large_data/                # Large binary files (Recommended for .gitignore)
    ├── beagle/                # Individual genotype likelihood matrices
    ├── mafs/                  # Point allele frequency tables
    ├── saf/                   # Population sample allele frequency distributions
    ├── thetas/                # Per-site nucleotide diversity estimates
    └── fst/                   # Binary FST indices
```

### Output Directories tracked by git

#### 1. `inputfiles/`
Contains formatted list files of BAM paths passed to ANGSD processes.
* **`bamlist.txt`**: Plain text list of all BAM files used across the full population dataset.
* **`${pop}.${era}.bamlist.txt`**: BAM file paths partitioned by specific population and sampling era.

#### 2. `sites/`
Contains filtered coordinate position files used to control site inclusion across downstream analyses.
* **`all_callable.pos`**: All callable genomic positions extracted from the input BED file, including monomorphic sites.
* **`all_snps.pos`**: Full set of polymorphic SNPs passing quality and $p$-value filters (`SNP_pval < 1e-6`).
* **`selected_loci.pos`**: Position coordinates of candidate loci identified under selection by ACER (`FDR <= cutoff`).
* **`callable_neutral.pos`** (`.bin`, `.idx`): All callable positions minus loci putatively under selection. Used by `ANGSD_DIVERSITY`.
* **`snps_neutral.pos`**: Polymorphic SNPs minus candidate selected loci.
* **`snps_pruned.pos`**: LD-pruned polymorphic SNP coordinates.
* **`snps_pruned_neutral.pos`**: Polymorphic SNPs that are both LD-pruned and neutral. Passed to Beagle subsetting for PCA/Admixture.

#### 3. `selection/`
Contains output tables and plots from the ACER iterative selection scan.
* **`iteration_summary.tsv`**: Round-by-round summary of effective population size ($N_e$) estimation and selection iterations.
* **`final_test_results.tsv`**: Per-SNP selection scan statistics, including raw $\chi^2$ test $p$-values and FDR-adjusted values.
* **`ne_bootstrap.tsv`**: Bootstrap replicates used to estimate baseline $N_e$.
* **`*.png`**: Manhattan plots and diagnostic visualizations generated during the selection scan.

#### 4. `pcangsd/`
Contains population structure and individual admixture results derived from LD-pruned neutral loci.
* **`${params.species}.cov`**: Estimated individual covariance matrix for PCA calculation.
* **`${params.species}.*.Q`**: Individual ancestry proportion matrix (Q-matrix) generated by PCAngsd admixture analysis.
* **`${params.species}.pcangsd.pdf`**: Publication-ready PCA scatter plot colored by population and era.
* **`${params.species}.admixture.pdf`**: Stacked bar plot displaying individual admixture proportions.

#### 5. `diversity/`
Summary nucleotide diversity files calculated across all callable neutral loci (monomorphic + polymorphic).
* **`${pop}_${era}.sfs`**: Folded site frequency spectrum (SFS) calculated across callable neutral sites.
* **`${pop}_${era}.pestPG`**: Windowed summary statistics generated by `thetaStat`, including pairwise nucleotide diversity ($\pi$), Watterson's $\theta$, and Tajima's $D$.

#### 6. `fst/`
FST between historical and modern populations derived from LD-pruned neutral loci.
* **`${pop}_hist_vs_mod.global_fst.txt`**: Global weighted and unweighted FST.
* **`${pop}_hist_vs_mod.2dsfs`**: 2D Site Frequency Spectrum matrix.


### `large_data/` Directory (Untracked / Ignored by git)

To exclude heavy intermediate data from Git, add `results/large_data/` to your `.gitignore` file.

| Subdirectory | File Pattern | Description |
| :--- | :--- | :--- |
| **`large_data/beagle/`** | `${params.species}.beagle.gz` | Full-genome individual genotype likelihood matrices in Beagle format. |
| **`large_data/mafs/`** | `${params.species}.mafs.gz`<br>`${pop}_${era}.mafs.gz` | Allele frequency tables for the entire dataset and individual population/era pairs. |
| **`large_data/saf/`** | `${pop}_${era}.saf.gz`<br>`${pop}_${era}.saf.idx`<br>`${pop}_${era}.saf.pos.gz` | Binary Sample Allele Frequency distributions used by `realSFS` to estimate diversity. |
| **`large_data/thetas/`** | `${pop}_${era}.thetas.gz`<br>`${pop}_${era}.thetas.idx` | Per-site binary likelihoods for nucleotide diversity estimates. |
| **`large_data/fst/`** | `${pop}_hist_vs_mod.fst.gz`<br>`${pop}_hist_vs_mod.fst.idx` | Binary FST indices. |