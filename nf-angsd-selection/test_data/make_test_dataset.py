#!/usr/bin/env python3
"""Generate a small synthetic ground-truth test dataset for nf-angsd-selection."""

import csv
import os
import random
import subprocess

random.seed(42)
OUT_DIR = "test_data"
os.makedirs(f"{OUT_DIR}/bams", exist_ok=True)

CONTIGS = {"chr1": 120_000, "chr2": 120_000}
REGIONS = ["RegionA", "RegionB"]
ERAS = ["historic", "modern"]
N_IND = 5  # 5 per region-era = 20 BAMs total
READ_LEN = 100
BASE_READS_PER_CONTIG = 7_200  # ~6x depth per sample

TRANSVERSIONS = {"A": "C", "C": "A", "G": "T", "T": "G"}
TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}

# 1. Generate reference FASTA + .fai
ref_seqs = {
    ctg: [random.choice("ACGT") for _ in range(length)]
    for ctg, length in CONTIGS.items()
}
ref_fa = f"{OUT_DIR}/ref.fa"
with open(ref_fa, "w") as f:
    for ctg, seq in ref_seqs.items():
        f.write(f">{ctg}\n")
        for i in range(0, len(seq), 80):
            f.write("".join(seq[i : i + 80]) + "\n")
subprocess.run(["samtools", "faidx", ref_fa], check=True)

# 2. Generate 1-based ANGSD BED file (masks chr1:100001-110000 as a repeat region)
bed_file = f"{OUT_DIR}/ref.repma.angsd.txt"
with open(bed_file, "w") as f:
    f.write("chr1\t1\t100000\n")
    f.write("chr1\t110001\t120000\n")
    f.write("chr2\t1\t120000\n")

with open(f"{OUT_DIR}/contig_list.txt", "w") as f:
    for ctg in CONTIGS:
        f.write(f"{ctg}\n")

# 3. Define ground-truth variants (0-based positions)
# Map: (ctg, pos0) -> dict of {(region, era): alt_freq} and metadata
variants = {}

def add_site(ctg, pos0, alt_base, freqs, tag):
    variants[(ctg, pos0)] = {
        "ref": ref_seqs[ctg][pos0],
        "alt": alt_base,
        "freqs": freqs,
        "tag": tag,
    }

# A) 160 Neutral transversion SNPs (80 per contig, outside masked/QC-fail windows)
import math
DRIFT_SD_SCALE = math.sqrt(1.0 - math.exp(-50.0 / (2.0 * 150.0)))  # Ne = 150, t = 50 gen
neutral_positions = {
    "chr1": [p for p in range(1_000, 99_000, 1_100) if not (15_000 <= p <= 16_000 or 20_000 <= p <= 21_500 or 30_000 <= p <= 31_500 or 45_000 <= p <= 46_000 or 75_000 <= p <= 76_000)][:80],
    "chr2": [p for p in range(1_000, 79_000, 900) if not (15_000 <= p <= 16_000 or 45_000 <= p <= 46_000 or 75_000 <= p <= 76_000)][:80],
}
for ctg, positions in neutral_positions.items():
    for pos0 in positions:
        p_anc = random.uniform(0.20, 0.50)
        # Regional divergence + small temporal drift
        p_regA = min(0.85, max(0.05, p_anc + random.gauss(0.10, 0.04)))
        p_regB = min(0.85, max(0.05, p_anc - random.gauss(0.10, 0.04)))
        sd_A = math.sqrt(p_regA * (1.0 - p_regA)) * DRIFT_SD_SCALE
        sd_B = math.sqrt(p_regB * (1.0 - p_regB)) * DRIFT_SD_SCALE
        freqs = {
            ("RegionA", "historic"): p_regA,
            ("RegionA", "modern"): min(0.95, max(0.05, p_regA + random.gauss(0, sd_A))),
            ("RegionB", "historic"): p_regB,
            ("RegionB", "modern"): min(0.95, max(0.05, p_regB + random.gauss(0, sd_B))),
        }
        add_site(ctg, pos0, TRANSVERSIONS[ref_seqs[ctg][pos0]], freqs, "neutral")

# B) 6 Strong selection sweep SNPs (p: 0.05 -> 0.95)
for ctg in ["chr1", "chr2"]:
    for pos0 in [15_499, 45_499, 75_499]:
        freqs = {
            ("RegionA", "historic"): 0.05,
            ("RegionA", "modern"): 0.95,
            ("RegionB", "historic"): 0.05,
            ("RegionB", "modern"): 0.95,
        }
        add_site(ctg, pos0, TRANSVERSIONS[ref_seqs[ctg][pos0]], freqs, "selected")

# C) 10 Tightly linked transversion SNPs (perfect LD block on chr1)
ld_positions = [20_099 + i * 100 for i in range(10)]
for pos0 in ld_positions:
    freqs = {
        ("RegionA", "historic"): 0.40,
        ("RegionA", "modern"): 0.40,
        ("RegionB", "historic"): 0.40,
        ("RegionB", "modern"): 0.40,
    }
    add_site("chr1", pos0, TRANSVERSIONS[ref_seqs["chr1"][pos0]], freqs, "ld_block")

# D) 10 Transition SNPs (should be filtered out by -rmTrans 1 / -noTrans 1)
for i in range(10):
    pos0 = 30_099 + i * 100
    freqs = {r_e: 0.35 for r_e in [("RegionA", "historic"), ("RegionA", "modern"), ("RegionB", "historic"), ("RegionB", "modern")]}
    add_site("chr1", pos0, TRANSITIONS[ref_seqs["chr1"][pos0]], freqs, "transition")

# E) 5 Repeat-masked SNPs inside chr1:100001-110000 (should be filtered out by BED)
for i in range(5):
    pos0 = 101_000 + i * 1_000
    freqs = {r_e: 0.40 for r_e in [("RegionA", "historic"), ("RegionA", "modern"), ("RegionB", "historic"), ("RegionB", "modern")]}
    add_site("chr1", pos0, TRANSVERSIONS[ref_seqs["chr1"][pos0]], freqs, "masked")

# 4. Simulate diploid genotypes & sorted/indexed BAMs
samplesheet_rows = []
qual_str = "I" * READ_LEN  # Phred Q40

for region in REGIONS:
    for era in ERAS:
        for idx in range(1, N_IND + 1):
            sample_id = f"{region}_{era[:4]}_{idx:02d}"
            bam_path = os.path.abspath(f"{OUT_DIR}/bams/{sample_id}.bam")
            samplesheet_rows.append([sample_id, bam_path, f"{region}_{era[:4]}", era, region])

            # Sample diploid alleles (0=ref, 1=alt)
            geno = {}
            ld_hap = (1 if random.random() < 0.40 else 0, 1 if random.random() < 0.40 else 0)
            for (ctg, pos0), info in variants.items():
                if info["tag"] == "ld_block":
                    geno[(ctg, pos0)] = ld_hap
                else:
                    p = info["freqs"][(region, era)]
                    geno[(ctg, pos0)] = (1 if random.random() < p else 0, 1 if random.random() < p else 0)

            # Build SAM stream and pipe into samtools sort
            sort_proc = subprocess.Popen(
                ["samtools", "sort", "-o", bam_path, "-"],
                stdin=subprocess.PIPE,
                text=True,
            )
            sort_proc.stdin.write("@HD\tVN:1.6\tSO:coordinate\n")
            for ctg, length in CONTIGS.items():
                sort_proc.stdin.write(f"@SQ\tSN:{ctg}\tLN:{length}\n")

            read_num = 0
            for ctg, length in CONTIGS.items():
                # Standard uniform reads
                for _ in range(BASE_READS_PER_CONTIG):
                    start0 = random.randint(0, length - READ_LEN)
                    # Plant low-coverage dropout in chr2:80000-82000 for 4/5 RegionA historic samples
                    if ctg == "chr2" and region == "RegionA" and era == "historic" and idx > 1:
                        if 79_900 <= start0 <= 82_000:
                            continue
                    read_num += 1
                    seq = list(ref_seqs[ctg][start0 : start0 + READ_LEN])
                    for offset in range(READ_LEN):
                        pos0 = start0 + offset
                        if (ctg, pos0) in geno:
                            a1, a2 = geno[(ctg, pos0)]
                            # Plant reference bias (60% ref) in historic heterozygotes
                            if a1 != a2:
                                ref_prob = 0.60 if era == "historic" else 0.50
                                allele = 0 if random.random() < ref_prob else 1
                            else:
                                allele = a1
                            if allele == 1:
                                seq[offset] = variants[(ctg, pos0)]["alt"]
                    flag = 0 if random.random() < 0.5 else 16
                    sort_proc.stdin.write(
                        f"{sample_id}_{read_num}\t{flag}\t{ctg}\t{start0 + 1}\t60\t{READ_LEN}M\t*\t0\t0\t{''.join(seq)}\t{qual_str}\tXO:i:1\tX0:i:1\tNM:i:0\n"
                    )

                # Plant high-depth paralog spike at chr2:90001-90500 (~100x extra depth)
                if ctg == "chr2":
                    for _ in range(500):
                        start0 = random.randint(90_000, 90_500 - READ_LEN)
                        read_num += 1
                        seq = "".join(ref_seqs[ctg][start0 : start0 + READ_LEN])
                        sort_proc.stdin.write(
                            f"{sample_id}_spike_{read_num}\t0\t{ctg}\t{start0 + 1}\t60\t{READ_LEN}M\t*\t0\t0\t{seq}\t{qual_str}\tXO:i:1\tX0:i:1\tNM:i:0\n"
                        )

            sort_proc.stdin.close()
            sort_proc.wait()
            subprocess.run(["samtools", "index", bam_path], check=True)

with open(f"{OUT_DIR}/samplesheet.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample", "bam", "pop", "era", "region"])
    w.writerows(samplesheet_rows)

# 5. Write lightweight local test Nextflow config
with open(f"{OUT_DIR}/test.config", "w") as f:
    f.write(f"""
params {{
    samplesheet   = "{os.path.abspath(OUT_DIR)}/samplesheet.csv"
    contigs       = "{os.path.abspath(OUT_DIR)}/contig_list.txt"
    reference     = "{os.path.abspath(OUT_DIR)}/ref.fa"
    bed_file      = "{os.path.abspath(OUT_DIR)}/ref.repma.angsd.txt"
    outdir        = "{os.path.abspath(OUT_DIR)}/results"
    species       = "Tst"
    generations   = 50
    n_boot        = 100
}}
""")

print(f"Created test dataset in ./{OUT_DIR}/")