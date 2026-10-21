#!/usr/bin/env bash
# =============================================================================
# 04_annotate.sh
# Part 5: annotate PASS variants with Ensembl VEP 116 (offline cache) and
#         build three prioritization tiers:
#           Tier 1  ClinVar pathogenic / likely pathogenic, MAX_AF < 1%
#           Tier 2  HIGH impact, MAX_AF < 1%, MANE Select transcript
#           Tier 3  SIFT deleterious + PolyPhen probably damaging, MAX_AF < 1%,
#                   allele balance >= 0.25
#
# Run from the project folder:   bash scripts/04_annotate.sh
# Conda environment:             gvp2-vep   (conda activate gvp2-vep first)
# Input:                         results/NA12878.pass.vcf.gz (from 03_call_filter.sh)
# =============================================================================

set -euo pipefail

# ---- 1. VEP cache, release 116, GRCh38 (one-time download, very large) -------
mkdir -p ~/.vep
wget -c -P ~/.vep https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
tar -xzf ~/.vep/homo_sapiens_vep_116_GRCh38.tar.gz -C ~/.vep

# ---- 2. Annotate; one consequence per variant, exon rank before APPRIS -------
vep --input_file results/NA12878.pass.vcf.gz \
  --output_file results/NA12878.vep.rankpick.vcf.gz --vcf --compress_output bgzip \
  --cache --offline --dir_cache ~/.vep --assembly GRCh38 \
  --fasta reference/GRCh38_no_alt.fa --fork 8 \
  --symbol --canonical --mane --biotype --hgvs --sift b --polyphen b \
  --check_existing --af_gnomade --af_gnomadg --max_af \
  --pick --pick_order mane_select,mane_plus_clinical,canonical,rank,appris,tsl,biotype,ccds,length \
  --stats_file logs/genome_vep_rankpick_summary.html 2> logs/genome_vep_rankpick.log

# ---- 3. Tier 1: ClinVar pathogenic / likely pathogenic, rare -----------------
bcftools +split-vep -i 'CLIN_SIG~"pathogenic"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%HGVSp\t%MAX_AF\t%CLIN_SIG[\t%GT]\n' \
  results/NA12878.vep.rankpick.vcf.gz 2>/dev/null > results/tier1_clinvar_rankpick.tsv
awk -F'\t' '$8 != "." && $8 < 0.01 && $9 ~ /(^|&)(likely_)?pathogenic/' \
  results/tier1_clinvar_rankpick.tsv > results/tier1_final_rankpick.tsv

# ---- 4. Tier 2: HIGH impact, rare, MANE Select transcript ---------------------
bcftools +split-vep -i 'IMPACT="HIGH"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%HGVSc\t%HGVSp\t%MAX_AF\t%CLIN_SIG\t%MANE_SELECT[\t%GT]\n' \
  results/NA12878.vep.rankpick.vcf.gz 2>/dev/null > results/tier2_high_all.tsv
awk -F'\t' '$9 == "." || $9 < 0.01' results/tier2_high_all.tsv > results/tier2_final.tsv
awk -F'\t' '$11 != "."' results/tier2_final.tsv > results/tier2_final_clean.tsv

# ---- 5. Tier 3: predicted damaging missense, rare -----------------------------
bcftools +split-vep -i 'IMPACT="MODERATE"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%HGVSp\t%SIFT\t%PolyPhen\t%MAX_AF\t%CLIN_SIG\t%MANE_SELECT[\t%GT]\n' \
  results/NA12878.vep.rankpick.vcf.gz 2>/dev/null > results/tier3_moderate_all.tsv
awk -F'\t' '$8 ~ /^deleterious\(/ && $9 ~ /^probably_damaging/ && ($10 == "." || $10 < 0.01)' \
  results/tier3_moderate_all.tsv > results/tier3_final.tsv

# ---- 6. Allele balance (VAF) check; remove Tier 3 calls with VAF < 0.25 -------
cut -f 1,2 results/tier1_final_rankpick.tsv results/tier2_final.tsv results/tier3_final.tsv \
  | sort -k1,1V -k2,2n | uniq > results/tier_positions.txt
bcftools query -T results/tier_positions.txt -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' \
  results/NA12878.vep.rankpick.vcf.gz \
  | awk -F'\t' '{split($6, a, ","); printf "%s\t%.2f\n", $0, a[2]/(a[1]+a[2])}' \
  > results/tier_allele_balance.tsv
awk -F'\t' '$7 < 0.25 {print $1"\t"$2"\t"}' results/tier_allele_balance.tsv > results/tier_lowvaf_positions.txt
grep -v -F -f results/tier_lowvaf_positions.txt results/tier3_final.tsv > results/tier3_final_clean.tsv

# ---- 7. Tier counts ---------------------------------------------------------
wc -l results/tier1_final_rankpick.tsv results/tier2_final_clean.tsv results/tier3_final_clean.tsv \
  | tee logs/part5_tier_counts.txt

echo "Annotation finished."
