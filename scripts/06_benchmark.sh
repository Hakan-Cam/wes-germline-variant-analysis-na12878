#!/usr/bin/env bash
# =============================================================================
# 06_benchmark.sh
# Part 7: benchmark the calls against GIAB HG001 v4.2.1 with hap.py,
#         inside the evaluation region (callable AND high-confidence):
#           run 1  original filter  QUAL<30 || DP<10 || MQ<50
#           run 2  tuned filter     QUAL<30 || DP<10 || MQ<40
#
# Run from the project folder:   bash scripts/06_benchmark.sh
# Conda environment:             gvp2-happy   (conda activate gvp2-happy first)
# Inputs:                        files prepared by 05_eval_regions.sh
# =============================================================================

set -euo pipefail

GIAB_VCF=truth/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
GIAB_BED=truth/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
EVAL=results/benchmark/eval_regions.bed
REF=reference/GRCh38_no_alt.fa

# ---- 1. hap.py, original filter (MQ<50) -------------------------------------------
hap.py ${GIAB_VCF} results/NA12878.norm.vcf.gz \
  -f ${GIAB_BED} -T ${EVAL} -r ${REF} \
  -o results/benchmark/happy_NA12878 --threads 4 2> logs/happy_mq50.log

# ---- 2. hap.py, tuned filter (MQ<40) ----------------------------------------------
hap.py ${GIAB_VCF} results/NA12878.norm.mq40.vcf.gz \
  -f ${GIAB_BED} -T ${EVAL} -r ${REF} \
  -o results/benchmark/happy_NA12878_mq40 --threads 4 2> logs/happy_mq40.log

# ---- 3. Side-by-side summary: recall, precision, F1 ---------------------------------
for run in happy_NA12878 happy_NA12878_mq40; do
  echo "== ${run}"
  cut -d, -f1-7,11,12,14 results/benchmark/${run}.summary.csv | column -s, -t
  echo
done | tee logs/happy_summaries.txt

echo "Benchmark finished."
