#!/usr/bin/env bash
# =============================================================================
# 02_qc_align.sh
# Part 2: read quality control (FastQC)
# Part 3: index the reference, align reads (bwa mem), sort, fix mate
#         information, mark duplicates, index, and collect alignment stats
#
# Run from the project folder:   bash scripts/02_qc_align.sh
# Conda environment:             gvp2
# Run time:                      several hours (bwa index ~40 min, bwa mem longest)
# Disk space needed:             ~50 GB while running (intermediates are removed)
# =============================================================================

set -euo pipefail

# ---- 1. Read quality control -------------------------------------------------
fastqc -t 2 --outdir qc reads/SRR098401_1.fastq.gz reads/SRR098401_2.fastq.gz 2>&1 | tee logs/fastqc.log

# ---- 2. Reference indexes ------------------------------------------------------
samtools faidx reference/GRCh38_no_alt.fa
bwa index reference/GRCh38_no_alt.fa 2>&1 | tee logs/bwa_index.log

# ---- 3. Align and coordinate-sort ---------------------------------------------
# Read group: ID = run, SM = sample, PL = platform, LB = library
bwa mem -t 24 -R '@RG\tID:SRR098401\tSM:NA12878\tPL:ILLUMINA\tLB:lib1' \
  reference/GRCh38_no_alt.fa \
  reads/SRR098401_1.fastq.gz reads/SRR098401_2.fastq.gz \
  2> logs/bwa_full.log \
  | samtools sort -@ 4 -m 2G -o results/NA12878.sorted.bam - 2> logs/sort_full.log

samtools index -@ 8 results/NA12878.sorted.bam
samtools flagstat -@ 8 results/NA12878.sorted.bam | tee logs/flagstat_full.txt

# ---- 4. Mark duplicates (collate -> fixmate -> sort -> markdup) --------------
samtools collate -@ 8 -o results/NA12878.collated.bam results/NA12878.sorted.bam
samtools fixmate -@ 8 -m results/NA12878.collated.bam results/NA12878.fixmate.bam
samtools sort -@ 8 -m 2G -o results/NA12878.fixmate.sorted.bam results/NA12878.fixmate.bam
samtools markdup -@ 8 -f logs/markdup_stats.txt results/NA12878.fixmate.sorted.bam results/NA12878.markdup.bam

samtools index -@ 8 results/NA12878.markdup.bam
samtools flagstat -@ 8 results/NA12878.markdup.bam | tee logs/flagstat_markdup.txt

# ---- 5. Remove large intermediate BAMs (~16 GB each) --------------------------
rm results/NA12878.sorted.bam results/NA12878.sorted.bam.bai \
   results/NA12878.collated.bam results/NA12878.fixmate.bam \
   results/NA12878.fixmate.sorted.bam

echo "Alignment finished: results/NA12878.markdup.bam"
