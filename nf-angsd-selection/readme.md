## Overview of nf-angsd-selection
This is in beta testing. It runs fine, though it may have undiscovered bugs!

This pipeline is for genotyping and basic pop gen analyses from temporal genomic lcWGS data at two time points. 

The pipeline
* calls probabilistic genotypes with ANGSD (excluding transitions)
* only considers sites that meet minimum individual and maximum depth limits in each of the population-eras (set as a fraction and a multiplier, respectively, of the full number of individuals and the total average depth across all individuals)
* runs (optional) chi-squared (and cmh if >1 region) selection scans from data across two time-points
* does (optional) ld-pruning
* calculates diversity estimates within the ANGSD framework on both all callable sites (`all_callable`) and putatively neutral callable sites (`callable_neutral`, including monomorphic sites), using any chi-sq or cmh < fdr-cutoff as evidence of selection 
* makes PCA and admixture plots on putatively neutral ld-pruned SNPs (`snps_pruned_neutral`)
* calculates FST and plots a sliding window FST and SFS curves for all callable sites (`all_callable`) and putatively neutral callable sites (`callable_neutral`).

It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster (use `-profile wahab`). The pipeline starts from BAM files, such as those output by nf-trim-generode. It was developed from nf-angsd-diversity by Malin Pinsky (August 2026) with plenty of help from Gemini.

**Important for WAHAB HPC:** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Installation

The pipeline can be installed along with the others in this repo directly from git

```bash

git clone https://github.com/mariannedehasque/nf-pipelines.git

```

## Required files

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* ANGSD-format BED file with masked repeat regions (eg, output by nf-trim-generode or by Generode)
* BAM files

In addition, we need to create the following files too (see instructions below):
* Samplesheet
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
* `pop`: currently unused by the workflow (retained for metadata compatibility).
* `era`: case-insensitive label indicating timepoint. Accepted values are `historic` (aliases: `historical`, `hist`, `era1`) and `modern` (aliases: `mod`, `era2`). Do not mix different aliases within the same `region` (e.g., use `historic` consistently rather than mixing `historic` and `historical`), and do not use arbitrary labels or sampling years—unrecognized labels are silently excluded from diversity comparisons, $F_{ST}$, and ACER.
* `region`: pairs historical and modern samples for direct temporal comparison. **Every `region` must contain samples from both eras**; regions with only one era are dropped from $F_{ST}$, diversity plots, and ACER selection scans.
* **Global site intersection:** Callable sites (`all_callable.pos`) are intersected across **all** `region` × `era` combinations in the samplesheet. Including a single low-coverage population/era will reduce the number of callable sites retained across the entire dataset.

If you want to use mapdamage-rescaled BAM files for historical samples (eg, from Generode or nf-trim-generode), point to those files in this input file.

Think carefully about which individuals to use. It may be useful to remove very low-depth individuals (e.g., <1x or <0.5x) based, for example, on `results/summary_reports/depth_by_individual_histogram.pdf` from `nf-trim-generode`. If low-depth individuals are included, consider reducing `min_ind_frac` so that lots of loci are not removed.


### BED file
The ANGSD-style BED file contains information on the sites that will be analyzed. We will use this file to remove repeat regions from the analysis. This files corresponds to the `-sites` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Sites).

This file is created automatically as part of the nf-trim-generode pipeline and output in `results/data/reference/*.repma.angsd.txt`.

Alternatively, to create the sites file from the GenErode output (after loading all modules):

```bash
awk '{print $1"\t"$2+1"\t"$3}' ./path/to/reference/reference.repma.bed > ./path/to/reference/reference.repma.angsd.txt

angsd sites index ./path/to/reference/reference.repma.angsd.txt
```

Think carefully about which contigs to use. It may be useful to remove very low-depth or high-depth contigs based, for example, on `results/summary_reports/depth_by_contig_histogram.pdf` from `nf-trim-generode`. These could be mtDNA, paralogs, or not from the target species.


### Contig file
Specify the contigs/scaffolds to include in the analysis, one per line (passed via `-r` to ANGSD in `QC_POP_SITES_CONTIG` and `ANGSD_GL_ALL`, and parallelized per contig in `LD_PRUNE_CONTIG`). Note that `ANGSD_GL_POP`, `ANGSD_DIVERSITY`, and `ANGSD_FST` run across all intersected sites jointly rather than submitting separate jobs per contig.

This code will pull a contig list out of the reference genome:
```
awk '/^>/ {sub(/^>/, ""); print $1}' /path/to/reference.fna > inputfiles/contig_list.txt
```

## Configuration

In the `main.nf` file, adjust the file paths for your project's reference and bed file.

In addition, the following parameters should be defined:
* `min_ind_frac`: The fraction of individuals that have to be represented at a locus/site to use it. This is applied per region-era and the intersection of the passing loci in each region-era is used (default 0.7)
* `max_depth_mult`: The multiple of the average depth to use as a ceiling in order to include a locus/site. This filter is applied per region-era and the intersection of the passing loci in each region-era is used (default 10)
* `ld_prune`: Whether to prune out SNPs in linkage disequilibrium within each contig (`ngsLD` + `prune_graph`, default `true`). Doing this is generally a good idea. When `false`, `snps_pruned.pos` defaults to `all_snps.pos` and PCA/Admixture run on unpruned neutral SNPs.
* `max_kb_dist`: LD pruning uses a sliding window to test for disequilibrium. This sets the window width. Narrower will run faster, wider will run more slowly (default 50)
* `min_weight`: Loci that are correlated (r2) more than this threshold will be trimmed out (default 0.2)
* `run_selection`: Whether to run the iterative ACER selection scan (default `true`). When `false`, `selected_loci.pos` is empty and `callable_neutral` / `snps_pruned_neutral` are identical to `all_callable` / `snps_pruned`.
* `generations`: Number of generations to use for Ne calculations (default 114 for Albatross to contemporary populations with a 1 year generation time)
* `fdr_cutoff`: False detection rate (FDR) cutoff for iterative selection scan (default 0.05)
* `max_rounds`: Maximum number of rounds to run for iterative selection scan (default 20)
* `n_boot`: Number of bootstrap replicates to run for Ne calculations (default 1000)
* `min_ind`: Minimum number of individuals to include in Ne calculations (default 4)
* `fst_window`: Window size for sliding window FST calculations (default 50000)
* `fst_step`: Step size for sliding window FST calculations (default 10000)


```bash
// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // sample sheet csv with columns for sample name, bam path, pop, era, and region. era is historic or modern. region matches historic and modern populations for diversity calculations. pop is not used
params.contigs     = "${projectDir}/inputfiles/contig_list.txt" // text file with one contig name per line. Contig names must match those in the reference genome FASTA file.
params.outdir      = "${projectDir}/results" // output directory for all results
params.reference   = "${projectDir}/data/reference/reference.fasta"
params.bed_file    = "${projectDir}/data/reference/reference.repma.angsd.txt"
params.species     = "Cvi" // 3-letter species code for output file naming. used in ANGSD output files and ACER results.
params.min_ind_frac   = 0.70  // Minimum percentage of individuals required per site (default 0.70)
params.max_depth_mult = 10    // Max depth threshold as a multiplier of the average total expected coverage (default 10)
params.ld_prune    = true   // Set to true to enable LD pruning (default true)
params.max_kb_dist = 50     // Maximum pairwise distance in kb to test for LD if pruning
params.min_weight  = 0.2    // Minimum r2 threshold for pruning filter (default 0.2)
params.run_selection = true  // Set to false to skip ACER by default (default true)
params.generations = 114    // Number of generations to use for Ne calculations. Default is 114 for Albatross to contemporary populations with a 1 year generation time. Adjust for other species as needed.
params.fdr_cutoff  = 0.05   // FDR cutoff for iterative selection scan (default 0.05)
params.max_rounds  = 20     // Maximum number of rounds to run for iterative selection scan (default 20)
params.n_boot      = 1000   // Number of bootstrap replicates to run for Ne calculations (default 1000)
params.min_ind     = 4      // Minimum number of individuals to include in Ne calculations (default 4)
params.fst_window = 50000  // Window size in bp for sliding window FST calculations (default 50000)
params.fst_step   = 10000  // Step size in bp for sliding window FST calculations (default 10000)
```

## Methodological & Filtering Assumptions

* **BAM indexing & depth threshold estimation:** All input BAMs must have `.bai`/`.csi` index files (`samtools idxstats` is used to count mapped reads). Expected mean depth is estimated from mapped read counts multiplied by the mean read length of the **first 1,000 reads** in each BAM, divided by total target bases in `bed_file` (`sum(end - start + 1)`, assuming non-overlapping 1-based inclusive intervals). `minInd` and `setMaxDepth` are rounded up to the nearest integer.
* **Hardcoded ANGSD & statistical filters:**
  * Base and mapping quality: `-minMapQ 25 -minQ 30 -uniqueOnly 1 -remove_bads 1`
  * Transitions excluded during genotype likelihood/frequency estimation (`-rmTrans 1` in `ANGSD_GL_ALL`; `-noTrans 1` in `ANGSD_GL_POP`).
  * Global SNP calling (`ANGSD_GL_ALL`): `-SNP_pval 1e-6 -skipTriallelic 0`, also subject to global `min_ind_frac` and `max_depth_mult` across the combined dataset.
  * Ancestral state & SFS folding: `ANGSD_GL_POP` supplies `-anc <reference.fasta>`, polarizing alleles against the reference assembly. `ANGSD_DIVERSITY` computes a **folded** 1D SFS (`-fold 1`) and summarizes $\pi$ (`tP / nSites`) across contigs/scaffolds with 10,000 bootstrap replicates (`PLOT_PI_DIVERSITY`). `ANGSD_FST` computes an **unfolded** 2D SFS (`realSFS` without `-fold 1`).
* **Neutral site masking (`selected_loci.pos`):** In `PREPARE_SITE_SETS`, a SNP is flagged as under selection if **either** the multi-region `cmh_fdr < fdr_cutoff` **or** any individual region's `<region>_chisq_fdr < fdr_cutoff`. Any locus flagged in at least one region is excluded from `callable_neutral` and `snps_pruned_neutral` across all regions.
* **PCAngsd $K$ selection:** `pcangsd --admix` automatically infers the number of ancestral populations ($K$) from the significant eigenvalues of the covariance matrix; $K$ is not manually fixed.

## To Run

It is nice to run nextflow pipelines in tmux so that you can keep it running in the background even when logging of.

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
# Load the nextflow module (version 23)
module load container_env nextflow/23.10.1

# Activate bash and your conda environment
bash

# Start the run
nextflow run main.nf -profile wahab -resume
```

The `nextflow.config` file holds the profiles you can choose, so that you can use different profiles to match different compute environments. At the moment, it only has a profile for `wahab`.

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

### Test dataset
A script to create a test dataset is included. It needs `samtools` available, so to run it on an interactive node on Wahab:
```
salloc # start the interactive node
bash
module load container_env samtools
mkdir -p /tmp/$USER/bin
echo -e '#!/bin/bash\nexec crun.samtools samtools "$@"' > /tmp/$USER/bin/samtools
chmod +x /tmp/$USER/bin/samtools
export PATH="/tmp/$USER/bin:$PATH"
python3 test_data/make_test_dataset.py
exit # to leave the salloc
```
This will create test data files in `test_data/`. The reference genome is 240,000 bp and has a 10,000 bp repeatitive region (chr1:100001–110000)plus a paralog on chr2. There are 191 SNPs, including 181 transversions (160 neutral + 6 sweeps + 10 in an LD-block + 5 in the repititive region)  and 10 transitions. Simulated Ne is 150. There is a 60% reference-allele bias planted in historic heterozygotes. Region A has higher ALT alleles frequencies and diversity than Region B. Both regions experience drift and lose some diversity.

To analyze the test dataset in a tmux shell:
```
tmux new -t nextflow
bash
module load container_env nextflow/23.10.1
nextflow run main.nf -profile wahab -c test_data/test.config -resume
```

Then verify the results:

1. Check site counts across all 7 libraries. It might not detect all six selected loci:
```
#    all_callable ~ 227,500 | all_snps == 176 | selected_loci == 6 | snps_pruned ~ 167
cat test_data/results/sites/site_counts.tsv
```

2. Verify only the 6 planted sweep SNPs (chr1/chr2: 15500, 45500, 75500) were detected. It might not detect all six:
```
cat test_data/results/large_data/sites/selected_loci.pos
```

3. Removing the 4 selected outliers after Round 1 should increase estimated $N_e$ and get closer to the simulated $N_e = 150$:
```
cat test_data/results/selection/iteration_summary.tsv
```

4. Compare global FST (all_callable > callable_neutral due to the 6 sweep SNPs):
```
head test_data/results/fst/*global_fst.txt
```

5. Verify historical reference bias is higher than modern by a few percentage points. We don't expect any strand bias (< 5%)
```
cat test_data/results/qc_bias/qc_bias_summary.tsv
```

6. Check that the [PCA](test_data/results/PCAngsd/Tst.pcangsd.plot.png) separates the regions along PC1 and the time points along PC2.

7. Check that $\pi$ is about 0.0003, which includes monomorphic sites, is higher in Region A than B, and is higher in historical than modern (by a little):
```
cat test_data/results/diversity/pi_summary_*.tsv
```
The callable neutral plots for [Region A](test_data/results/diversity/pi_historic_vs_modern_RegionA_callable_neutral.png) and [Region B](test_data/results/diversity/pi_historic_vs_modern_RegionB_callable_neutral.png) show this as well.

## Output
## Pipeline Output Directory Structure

All pipeline outputs are written to the directory specified by `--outdir` (default: `./results/`). Large files are isolated in the `large_data/` directory so they can be easily managed or excluded via `.gitignore`.

```
results/
├── diversity/                 # Summary nucleotide diversity statistics and comparison plots
├── fst/                       # Global FST values, SFS plot, and Manhattan plots
├── inputfiles/                # Input BAM file lists
├── PCAngsd/                   # Population structure (PCA & Admixture)
├── qc_bias/                   # Sequencing bias diagnostics (reference bias, strand bias)
├── selection/                 # ACER iterative selection scan outputs
├── sites/                     # Summary metrics across generated site set libraries
└── large_data/                # Large binary, intermediate, and site library files (Ignored by git)
    ├── beagle/                # Individual genotype likelihood matrices
    ├── counts/                # Per-site, per-individual allele count matrices
    ├── fst/                   # Binary FST indices, 2D SFS matrices, and windowed FST tables
    ├── mafs/                  # Point allele frequency tables
    ├── saf/                   # Population sample allele frequency distributions
    ├── sites/                 # Base coordinate indices, region mappings, and processed site libraries
    ├── snpstat/               # Per-site SNP and strand-bias test statistics
    └── thetas/                # Per-site nucleotide diversity estimates
```

### Output Directories Tracked by Git

#### 1. `inputfiles/`
Contains formatted list files of BAM paths passed to ANGSD processes
* **`bamlist.txt`**: Plain text list of all BAM files used across the full population dataset.
* **`${pop}.${era}.bamlist.txt`**: BAM file paths partitioned by specific population/region and sampling era.

#### 2. `sites/`
Contains site library summary metrics.
* **`site_counts.tsv`**: Summary table detailing total site counts across all processed site subsets (`all_callable`, `all_snps`, `selected_loci`, `callable_neutral`, `snps_neutral`, `snps_pruned`, `snps_pruned_neutral`) generated by `PREPARE_SITE_SETS`.

#### 3. `selection/`
Contains output tables and plots from the ACER iterative selection scan.
* **`iteration_summary.tsv`**: Round-by-round summary of effective population size ($N_e$) estimation and selection iterations.
* **`final_test_results.tsv`**: Per-SNP selection scan statistics, including raw $\chi^2$ test $p$-values and FDR-adjusted values.
* **`ne_bootstrap.tsv`**: Bootstrap replicates used to estimate baseline $N_e$.
* **`*.png`**: Manhattan plots and diagnostic visualizations generated during the selection scan.

#### 4. `PCAngsd/`
Contains population structure and individual admixture results derived from LD-pruned neutral loci.
* **`${params.species}.cov`**: Estimated individual covariance matrix for PCA calculation.
* **`${params.species}.*.Q`**: Individual ancestry proportion matrix (Q-matrix) generated by PCAngsd admixture analysis.
* **`${params.species}.pcangsd.plot.png`**: Publication-ready PCA scatter plot colored by population and era.
* **`${params.species}.admixture.png`**: Stacked bar plot displaying individual admixture proportions.

#### 5. `diversity/`
Summary nucleotide diversity files calculated across all callable loci (`all_callable`) and putatively neutral callable loci (`callable_neutral`), including both monomorphic and polymorphic sites, indicated by `${label}`.
* **`${region}_${era}_${label}.sfs`**: Folded site frequency spectrum (SFS) calculated across the specified callable site set.
* **`${region}_${era}_${label}.pestPG`**: Windowed summary statistics generated by `thetaStat`, including pairwise nucleotide diversity ($\pi$), Watterson's $\theta$, and Tajima's $D$.
* **`pi_summary_${region}_${label}.tsv`**: Summary table of mean nucleotide diversity ($\pi / \theta_P$) with 95% bootstrap confidence intervals across eras.
* **`pi_historic_vs_modern_${region}_${label}.png`**: Comparison plot of historic vs. modern nucleotide diversity per region and site set.

#### 6. `fst/`
$F_{ST}$ and 2D SFS metrics between historical and modern populations derived from all callable loci (`all_callable`) and putatively neutral callable loci (`callable_neutral`), including monomorphic sites, indicated by `${label}`.
* **`${region}_${label}_hist_vs_mod.global_fst.txt`**: Global unweighted and weighted $F_{ST}$ values, in that order.
* **`${region}_${label}_fst_manhattan.png`**: Manhattan plot displaying genome-wide $F_{ST}$ sliding window calculations across chromosomes.
* **`sfs_summary_all_datasets.png`**: Site frequency spectrum (SFS) comparison plot of historical vs. modern allele frequencies across all regions and callable site subsets (`all_callable` and `callable_neutral`). Useful for identifying bioinformatic problems (paralogs can appear as an excess of alleles at 50% frequency; jagged SFS are also suspect) or a loss of low-frequency alleles through time indicative of a bottleneck.

#### 7. `qc_bias/`
Quality-control diagnostics comparing sequencing and mapping biases between historical and modern samples across called SNPs (`all_snps.pos`), generated by `EVALUATE_QC_BIAS`.
* **`qc_bias_summary.tsv`**: Summary table of mean reference allele balance (`Ref / (Ref + Alt)`), and percentage of loci with significant strand bias ($p < 0.05$) per region and era.
* **`qc_bias_plots.png`**: Multi-panel diagnostic plot comparing reference bias density distributions and strand bias across populations and eras.
---

### `large_data/` Directory (Untracked / Ignored by Git)

Contains heavy intermediate data, binary index files, and site set libraries (excluded by adding `results/large_data/` to `.gitignore`).

#### 1. `large_data/sites/`
Base genomic coordinate indices, region mappings, and processed site set libraries.
* **`all_sites.pos`** (`.bin`, `.idx`, `.regions`): Reference BED genomic targets (all 1-based genomic positions extracted directly from the input BED file prior to sequence read coverage or quality filtering) generated by `INDEX_BED_SITES`.
* **`${params.species}.sites.txt`** (`.bin`, `.idx`, `.regions`, `.chrs`): Polymorphic site indices, region mappings, and chromosome lists output by `ANGSD_EXTRACT_SITES`. Not filtered by minind or maxdepth.
* **`all_callable.pos`** (`.bin`, `.idx`): QC-filtered shared callable sites (the intersected subset of reference BED positions passing `minind_frac` and `maxdepth_mult` thresholds across every region-era subset) extracted and indexed by `PREPARE_SITE_SETS`. Includes monomorphic sites.
* **`all_snps.pos`**: Full set of polymorphic SNPs passing quality and $p$-value filters.
* **`selected_loci.pos`**: Position coordinates of candidate loci identified under selection by ACER.
* **`callable_neutral.pos`** (`.bin`, `.idx`): All callable positions minus candidate selected loci.
* **`snps_neutral.pos`**: Polymorphic SNPs minus candidate selected loci.
* **`snps_pruned.pos`**: LD-pruned polymorphic SNP coordinates.
* **`snps_pruned_neutral.pos`** (`.bin`, `.idx`): Polymorphic SNPs that are both LD-pruned and neutral.

#### 2. `large_data/beagle/`
Full-genome individual genotype likelihood matrices in Beagle format.
* **`${params.species}.beagle.gz`**: Concatenated individual genotype likelihood matrix generated by `ANGSD_COLLECT_OUTPUT`. The positions are those in `all_snps.pos`.

#### 3. `large_data/mafs/`
Point allele frequency tables across populations and subsets.
* **`${params.species}.mafs.gz`**: Combined SNPS frequency table across all individuals. The positions are those in `all_snps.pos`.
* **`${region}_${era}.mafs.gz`**: Major/minor allele frequency outputs per population/era subset, including monomorphic. The positions are those in `all_callable.pos`.

#### 4. `large_data/counts/`
* **`${region}_${era}.counts.gz`**: Per-individual A, C, G, T base counts across `all_callable.pos` (`-dumpCounts 4`), used by `EVALUATE_QC_BIAS` to compute reference allele balance.

#### 5. `large_data/snpstat/`
* **`${region}_${era}.snpStat.gz`**: ANGSD per-site quality and strand-bias test statistics (`-dosnpstat 1`) across `all_callable.pos`, used by `EVALUATE_QC_BIAS`.

#### 6. `large_data/saf/`
Population sample allele frequency distributions used for downstream site frequency spectra.
* **`${region}_${era}.saf.gz`** (`.idx`, `.pos.gz`): Binary Sample Allele Frequency distributions generated by `ANGSD_GL_POP`. The positions are those in `all_callable.pos`.

#### 7. `large_data/thetas/`
Per-site binary likelihoods for nucleotide diversity estimates.
* **`${region}_${era}_${label}.thetas.gz`** (`.idx`): Binary per-site theta likelihood output generated by `realSFS saf2theta`. The positions correspond to `all_callable.pos` and `callable_neutral.pos` (indicated by `${label}`).

#### 8. `large_data/fst/`
Binary $F_{ST}$ indices, 2D-SFS matrices, and sliding window tables.
* **`${region}_${label}_hist_vs_mod.fst.gz`** (`.idx`): Binary per-site $F_{ST}$ index files generated by `realSFS fst index`. The positions correspond to `all_callable.pos` and `callable_neutral.pos` (indicated by `${label}`).
* **`${region}_${label}_hist_vs_mod.2dsfs`**: 2D Site Frequency Spectrum matrix estimated between historical and modern populations. The positions correspond to `all_callable.pos` and `callable_neutral.pos` (indicated by `${label}`).
* **`${region}_${label}_windowed_fst.txt`**: Genome-wide sliding window $F_{ST}$ table output by `realSFS fst stats2`. The positions correspond to `all_callable.pos` and `callable_neutral.pos` (indicated by `${label}`).