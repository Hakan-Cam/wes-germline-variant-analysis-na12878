#!/usr/bin/env bash
# =============================================================================
# 05_eval_regions.sh
# Part 6: GIAB truth statistics, callable regions from the BAM (>=10x, MAPQ>=20),
#         evaluation region = callable AND GIAB high-confidence, sanity counts
# Part 7: which filter condition removed true variants; re-filter with MQ<40
#
# Run from the project folder:   bash scripts/05_eval_regions.sh
# Conda environment:             gvp2
# Inputs:                        results/NA12878.markdup.bam, results/NA12878.norm.vcf.gz,
#                                results/NA12878.pass.vcf.gz, GIAB files in truth/
# =============================================================================

set -euo pipefail

mkdir -p results/coverage results/benchmark

GIAB_VCF=truth/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
GIAB_BED=truth/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
EVAL=results/benchmark/eval_regions.bed

# ---- 1. Truth set statistics ---------------------------------------------------
bcftools stats ${GIAB_VCF} > logs/giab_truth_stats.txt

# ---- 2. Callable regions: depth >= 10x with MAPQ >= 20 ---------------------------
mosdepth --threads 4 --no-per-base --mapq 20 --quantize 0:1:10: \
  results/coverage/NA12878 results/NA12878.markdup.bam
zcat results/coverage/NA12878.quantized.bed.gz \
  | awk -v OFS='\t' '$4 == "10:inf" {print $1, $2, $3}' > results/coverage/NA12878.callable.bed

# ---- 3. Evaluation region: callable AND GIAB high-confidence -------------------
bedtools intersect -a results/coverage/NA12878.callable.bed -b ${GIAB_BED} > ${EVAL}

# ---- 4. Sanity counts inside the evaluation region -----------------------------
bcftools stats -T ${EVAL} ${GIAB_VCF} > logs/giab_truth_eval_stats.txt
bcftools stats -T ${EVAL} results/NA12878.pass.vcf.gz > logs/pass_eval_stats.txt

# ---- 5. Tabix index of the query VCF for hap.py ---------------------------------
bcftools index -f -t results/NA12878.norm.vcf.gz

# ---- 6. Filter diagnosis: true variants labelled LowQual by the MQ<50 filter -----
bcftools view -f LowQual -T ${EVAL} -Oz -o results/benchmark/lowqual_eval.vcf.gz results/NA12878.norm.vcf.gz
bcftools index -f -t results/benchmark/lowqual_eval.vcf.gz
bcftools isec -n=2 -w1 -Oz -o results/benchmark/lowqual_in_truth.vcf.gz \
  results/benchmark/lowqual_eval.vcf.gz ${GIAB_VCF}

LQ=results/benchmark/lowqual_in_truth.vcf.gz
{
  echo "True variants labelled LowQual:  $(bcftools view -H ${LQ} | wc -l)"
  echo "  failing QUAL<30:               $(bcftools view -H -i 'QUAL<30' ${LQ} | wc -l)"
  echo "  failing INFO/DP<10:            $(bcftools view -H -i 'INFO/DP<10' ${LQ} | wc -l)"
  echo "  failing INFO/MQ<50:            $(bcftools view -H -i 'INFO/MQ<50' ${LQ} | wc -l)"
  echo "  with MQ 40-50 (rescuable):     $(bcftools view -H -i 'INFO/MQ>=40 && INFO/MQ<50' ${LQ} | wc -l)"
} | tee logs/lowqual_in_truth_conditions.txt

# ---- 7. Re-filter with the GATK-recommended MQ<40 threshold ------------------------
bcftools filter -m x -s LowQual -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40' \
  -Oz -o results/NA12878.norm.mq40.vcf.gz results/NA12878.norm.vcf.gz
bcftools index -f -t results/NA12878.norm.mq40.vcf.gz

echo "Evaluation region and query files ready for hap.py."
