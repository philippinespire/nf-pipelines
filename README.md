## Overview
Repository to document custom nextflow pipelines. Pipelines were written for the WAHAB cluster.

Documentation and testing is work in progress.

## Content

```
nf-pipelines/
├── README.md
├── docs/
│   ├── lab_notebook.md
├── nf-angsd-diversity/         # Runs angsd for genotypes, calculates diversity and a PCA.
├── nf-trim-generode/           # Takes historical and modern reads. Trims modern to average historical length and maps both with the same mapper. Also masks repeats in the reference. A merger of Generode and nf-trim-merged-unmerged.
└── nf-trim-merged-unmerged/    # Takes modern reads. Trims modern to an average length and maps. Usually run after Generode.


```
