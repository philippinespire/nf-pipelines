#!/usr/bin/env bash
set -euo pipefail

# One-time setup: build BWA index files for the synthetic test reference genome.
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_FASTA="${TEST_DIR}/reference/test_ref.fasta"

if [[ -f "${REF_FASTA}.bwt" ]]; then
    echo "BWA index already exists for ${REF_FASTA}. Skipping."
    exit 0
fi

echo "Indexing ${REF_FASTA} with BWA..."
module load container_env bwa
crun bwa index "${REF_FASTA}"

echo "Done! Launch the test run with: nextflow run main.nf -profile test,wahab -resume"