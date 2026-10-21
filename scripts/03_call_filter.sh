#!/usr/bin/env bash
# =============================================================================
# 03_call_filter.sh
# Part 4: call variants per chromosome (bcftools mpileup | call), combine,
#         soft-filter, normalize, and collect summary statistics
# Part 5: keep PASS calls, excluding chrY and chrM (female sample)
#
# Run from the project folder:   bash scripts/03_call_filter.sh
# Conda environment:             gvp2
# Input:                         results/NA12878.markdup.bam (from 02_qc_align.sh)
# =============================================================================

set -euo pipefail

mkdir -p results/per_chr

# ---- 1. Call variants: one chromosome per background job ---------------------
for c in chr{1..22} chrX chrY chrM; do
  ( bcftools mpileup -Ou -f reference/GRCh38_no_alt.fa -r $c \
      -q 20 -Q 20 -d 1000 -a FORMAT/AD,FORMAT/DP results/NA12878.markdup.bam \
    | bcftools call -mv -Oz -o results/per_chr/NA12878.$c.vcf.gz ) 2> logs/call_$c.log &
done
wait

# ---- 2. Combine chromosomes in genome order -----------------------------------
bcftools concat -Oz -o results/NA12878.raw.vcf.gz \
  results/per_chr/NA12878.chr{1..22}.vcf.gz results/per_chr/NA12878.chr{X,Y,M}.vcf.gz
bcftools index results/NA12878.raw.vcf.gz

# ---- 3. Soft filter: label low-quality calls, keep everything -----------------
bcftools filter -s LowQual -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<50' \
  -Oz -o results/NA12878.filtered.vcf.gz results/NA12878.raw.vcf.gz

# ---- 4. Normalize: left-align indels, split multiallelic sites ----------------
bcftools norm -f reference/GRCh38_no_alt.fa -m -any \
  -Oz -o results/NA12878.norm.vcf.gz results/NA12878.filtered.vcf.gz
bcftools index results/NA12878.norm.vcf.gz

# ---- 5. Summary statistics of PASS calls --------------------------------------
bcftools stats -f PASS results/NA12878.norm.vcf.gz > logs/genome_stats_pass.txt
bcftools view -H -f PASS results/NA12878.norm.vcf.gz | cut -f 1 | sort -V | uniq -c | tee logs/pass_per_chrom.txt

# ---- 6. PASS set for annotation (Part 5); chrY/chrM removed -------------------
bcftools view -f PASS -t ^chrY,chrM -Oz -o results/NA12878.pass.vcf.gz results/NA12878.norm.vcf.gz

echo "Calling finished: results/NA12878.norm.vcf.gz and results/NA12878.pass.vcf.gz"
