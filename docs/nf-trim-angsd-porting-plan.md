# nf-trim-angsd Porting Plan

## Goal

Extend the raw-read workflow so it can serve as the preprocessing front end for ANGSD-based analyses while borrowing the historical-read handling pattern from GenErode.

The target state is:

- modern reads are trimmed and mapped as they are today,
- the modern trimming length is computed from the average length of the historical reads,
- historical reads are not trimmed,
- all historical reads can optionally be rescaled with mapDamage,
- both branches converge on analysis-ready BAMs that can feed the downstream ANGSD pipeline.

## Current Baseline

The current [nf-trim-merged-unmerged](../nf-trim-merged-unmerged/main.nf) workflow already does most of the raw-read work:

- fastp on merged and unmerged reads,
- length trimming with seqtk for unmerged reads,
- mapping with BWA,
- duplicate removal,
- BAM merging,
- indel realignment,
- BAM QC,
- optional mapDamage,
- optional AMBER prep/run.

GenErode adds the historical-aware behavior we want to borrow:

- reference-repeat preparation that can be run inside the workflow,
- a dedicated historical mapDamage branch,
- optional rescaling of historical BAMs for downstream analyses,
- explicit configuration for which samples are historical.

## Migration Principles

1. Keep the ANGSD analysis pipeline separate from preprocessing.
2. Treat raw-read preprocessing as a shared upstream workflow.
3. Make the modern and historical paths explicit in configuration, not implicit in sample naming.
4. Avoid trimming historical reads unless a future requirement proves that is necessary.
5. Make mapDamage rescaling an all-or-none option for the historical cohort.
6. Keep merged and unmerged read handling separate so the workflow preserves the current raw-read structure.

## Proposed Workflow Split

### Modern branch

- Input: modern FASTQ pairs.
- Steps:
  - trim / merge / split as currently done, using a trim length derived from the historical-read length distribution,
  - map to the reference,
  - remove duplicates,
  - merge BAMs per biological sample if needed,
  - indel realignment,
  - BAM QC,
  - optional downstream use in ANGSD.

### Historical branch

- Input: historical FASTQ pairs.
- Steps:
  - skip trimming entirely,
  - optionally run FastQC before mapping,
  - map reads directly,
  - remove duplicates,
  - merge BAMs per biological sample if needed,
  - indel realignment,
  - optional mapDamage run on the realigned BAMs,
  - emit both the original BAM and the rescaled BAM for downstream analyses when rescaling is enabled.

## Functionality To Port From GenErode

### 1. Reference and repeat-mask preparation

Port the reference-prep and repeat-identification logic first.

Why this matters:

- it standardizes the reference used by both mapping and site filtering,
- it creates the repeat-masked regions that ANGSD should exclude,
- it removes manual pre-run setup steps.

### 2. Historical sample routing

Add explicit sample metadata for at least:

- sample ID,
- era or sample class,
- whether the sample is historical,
- whether the sample should be used to define the historical read-length average.

This should allow the workflow to decide whether a read pair enters the modern trim branch or the historical no-trim branch.

The historical sample list should be the only sample-level selector needed for the historical mapDamage path. If historical rescaling is enabled, it should apply to the full historical cohort rather than a subset.

### 3. Historical read-length summary

Add a preprocessing step that measures read lengths in the historical FASTQ inputs and computes the average length.

That value should drive the modern trimming length so the modern branch aligns with the historical read geometry rather than a hard-coded constant.

### 4. Optional mapDamage for historical BAMs

Add a dedicated historical BAM branch that can run mapDamage after mapping and deduplication.

The key output should be twofold:

- a mapDamage report for damage assessment,
- a rescaled BAM for the historical cohort when rescaling is enabled.

This mirrors the GenErode behavior where mapDamage is not just diagnostic; it can also produce BAMs for downstream use.

When enabled, mapDamage rescaling should be applied to all historical samples, and both the original and rescaled BAMs should remain available for downstream use.

### 5. Coverage and QC checkpoints

Keep or strengthen the depth/QC stage so the pipeline can surface:

- mean depth per sample,
- coverage after duplicate removal and realignment,
- whether the historical branch changes coverage distribution after rescaling.

### 6. Optional subsampling

If coverage balancing becomes important for ANGSD or population comparisons, port the subsampling concept as an opt-in step after BAM cleanup and before downstream analysis.

## Functionality To Leave Out For Now

These are useful in GenErode, but they should not be part of the first port:

- genotyping and VCF-centric filtering,
- CpG identification,
- ROH,
- snpEff,
- GERP,
- autosome/sex-chromosome BED generation,
- any analysis track that assumes variant calling is the primary output.

## Implementation Tasks

### Task 1: Internalize reference-repeat preparation

Move reference prep and repeat-mask generation into the pipeline so the workflow no longer depends on manual pre-run setup.

Acceptance criteria:

- the reference preparation step runs from the pipeline inputs,
- repeat-masked regions are emitted as pipeline outputs,
- downstream steps can consume the generated mask directly.

### Task 2: Add sample metadata for modern vs historical routing

Define a sample sheet or config structure that marks each sample as modern or historical.

Acceptance criteria:

- the workflow can split inputs into modern and historical branches,
- the historical sample list is sufficient to identify the historical cohort,
- no separate per-sample rescaling list is needed.

### Task 3: Measure historical read length and derive the modern trim length

Add a preprocessing step that calculates the average read length from the historical FASTQ inputs and uses that value to set the modern trimming length.

Acceptance criteria:

- historical reads are summarized before modern trimming starts,
- the modern trim length is derived from the historical average,
- the chosen trim length is visible in logs or outputs for traceability.

### Task 4: Split modern and historical raw-read handling

Preserve merged and unmerged read handling as separate streams and route them through the appropriate branch.

Acceptance criteria:

- merged and unmerged reads remain distinct in the workflow graph,
- modern reads are trimmed before mapping,
- historical reads bypass trimming
- historical reads have merged and unmerged tracks, like modern reads

### Task 5: Decide on which mapper to use based on trim length

Choose BWA ALN is the modern trim length is <= 80 bp, and choose BWA MEM is trim length is larger. The nf-pipelines/nf-trim-merged-unmerged/main.nf script has an example for how to choose among these two options. Both modern and historical reads should be mapped with the same mapper and parameters.

Acceptance criteria:

- bwa aln is used if trim length is <= 80 bp, and bwa mem is used otherwise
- both modern and historical reads are mapped in the same way

### Task 6: Add optional historical FastQC

Make FastQC available as an optional pre-mapping step for historical reads.

Acceptance criteria:

- the option can be enabled or disabled by configuration,
- enabling it does not force trimming,
- FastQC outputs are published alongside other preprocessing reports.

### Task 7: Add historical mapDamage support

Run mapDamage on historical BAMs after mapping and duplicate removal, with optional rescaling for the full historical cohort.

Acceptance criteria:

- a historical mapDamage flag controls whether the step runs,
- when enabled, all historical samples are processed,
- both original and rescaled historical BAMs are emitted for downstream use,
- rescaling is all-or-none for the historical cohort.

### Task 8: Preserve QC and depth reporting

Keep or strengthen the QC and coverage reporting so the new branch structure remains observable.

Acceptance criteria:

- per-sample depth summaries remain available,
- BAM QC runs after cleanup/realignment,
- the outputs make it clear which BAM variant is being summarized.

### Task 9: Wire the cleaned BAMs into ANGSD

Define the exact BAM artifact that downstream ANGSD modules should consume for each sample class.

Acceptance criteria:

- modern samples feed ANGSD from the cleaned modern BAMs,
- historical samples can feed ANGSD from either original or rescaled BAMs according to configuration,
- the preprocessing pipeline remains separate from ANGSD analysis modules.

## Execution Order

1. Internalize reference-repeat preparation.
2. Add sample metadata for modern vs historical routing.
3. Measure historical read length and derive the modern trim length.
4. Split modern and historical raw-read handling while keeping merged and unmerged streams separate.
5. Choose mapper based on trim length
6. Add optional historical FastQC.
7. Add historical mapDamage support and emit both historical BAM variants.
8. Preserve QC and depth reporting.
9. Wire the cleaned BAMs into ANGSD.

## Definition Of Done

- Modern samples trim, map, deduplicate, realign, and QC using a trim length derived from the historical read-length average.
- Historical samples skip trimming, can optionally run FastQC, map, deduplicate, realign, and optionally undergo mapDamage rescaling.
- Reference-repeat preparation runs inside the pipeline.
- Merged and unmerged reads remain separate in the workflow.
- bwa mem is used if trim length is >80, and bwa aln is used otherwise
- Both original and rescaled historical BAMs are available when rescaling is enabled.
- ANGSD consumes the appropriate cleaned BAMs downstream.
