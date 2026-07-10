#!/bin/bash -le
#SBATCH --job-name=collect_stats
#SBATCH -o collect_stats-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --array=1-292%100

## information

# this script collects and calculates summary statistics from modern samples in the GenErode pipeline.
# script is based on the historical version, so some results may not make a lot of sense.
# the script is run per sample, and adds output to the file "modern_sequencing_stats.txt"
# in order to get a sample list ("samp_list")from the metadata table, you can use sed as below:
# awk 'NR>1 {print$1}' ../../../../../jselwyn/spratelloides_delicatulus_lcwgs/GenErode_20k/config/modern_samples.txt | cut -d "_" -f1 | sort | uniq > samp_list
# the script is run within the main GenErode directory
# stats based on mapping quality (MQ) 25. 

## edits to make before running:
# change the {ref} to match your reference genome name
# change the path to your reference bed file
# change the path to the scripts directory

## usage of the script
# sbatch collect_stats_modern.sh

SAMPLES_FILE=/archive/carpenterlab/pire/marianne/spratelloides_delicatulus_lcwgs/QC/scripts/inputfiles/samp_list_modern
SAMPLENAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES_FILE)

echo ${SAMPLENAME}

## directories and files

# Generode directory
generode='/archive/carpenterlab/pire/jselwyn/spratelloides_delicatulus_lcwgs/GenErode_20k'
# Directory for output stats
outdir='/archive/carpenterlab/pire/marianne/spratelloides_delicatulus_lcwgs/QC/results/generode_stats'
# reference genome name
ref=GCA_049181635.1.Sde20k
# directory containing raw sequencing stats
raw_reads='data/raw_reads_symlinks/modern/stats'
# directory containing mapped reads stats 
mapped="results/modern/mapping/${ref}"
# path to reference repeat masked bed file
bed="reference/${ref}.bed"
# path to directory containing calculate.awk
scripts='/archive/carpenterlab/pire/marianne/spratelloides_delicatulus_lcwgs/QC/scripts'

## load software

module load container_env
module load samtools

## Change to generode directory
cd $generode

## summary statistics

shopt -s nullglob
# Sum total sequences across all R1 fastqc reports for this sample
#total_seqs=$(awk '{s+=$3} END{print (s==""?0:s)}' \
#  $raw_reads/${SAMPLENAME}_*_R1_fastqc/fastqc_data.txt)

# Sum accross R1 and R2 since most will not merge
total_seqs=$(awk '{s+=$3} END{print (s==""?0:s)}' \
  $raw_reads/${SAMPLENAME}_*_R1_fastqc/fastqc_data.txt \
  $raw_reads/${SAMPLENAME}_*_R2_fastqc/fastqc_data.txt)

# Sum mapped reads across all sorted.bam stats for this sample
total_mapped=$(awk 'FNR==5 {s+=$1} END{print (s==""?0:s)}' \
  $mapped/stats/bams_sorted/${SAMPLENAME}_*.sorted.bam.stats.txt)

# Get unique reads from rmdup_bam stats

total_uniq=$(awk 'NR==5 {print$1}' $mapped/stats/bams_rmdup/${SAMPLENAME}_*.merged.rmdup.bam.stats.txt)

# get mq 25 reads

total_MQ25=$(crun.samtools samtools view -c -q 25 $mapped/$SAMPLENAME.merged.rmdup.merged.realn.bam)

# get endogenous

endogenous=$(awk "BEGIN {printf \"%.5f\", $total_mapped / $total_seqs}")

# get complexity

complexity=$(awk "BEGIN {printf \"%.5f\", $total_uniq / $total_mapped}")

# get genome recovery rate

grr=$(awk "BEGIN {printf \"%.5f\", $total_MQ25 / $total_seqs}")

# get total coverage

total_cov=$(crun.samtools samtools depth -a -q 30 -Q 25 -b $bed $mapped/$SAMPLENAME.merged.rmdup.merged.realn.bam | awk '{{sum+=$3}} END {{ print sum/NR }}' | awk '{print$1}')

# get read lengths
# prints min max median mean

read_length=$(crun.samtools samtools view -q 25 $mapped/$SAMPLENAME.merged.rmdup.merged.realn.bam | awk '{print length($10)}' | awk -f $scripts/calculate.awk)

# print file

echo $SAMPLENAME $total_seqs $total_mapped $total_uniq $total_MQ25 $endogenous $complexity $grr $total_cov $read_length >> $outdir/modern_sequencing_stats.txt
