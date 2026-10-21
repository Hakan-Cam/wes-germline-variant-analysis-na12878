#!/usr/bin/env bash
# =============================================================================
# 01_download.sh
# Part 2: download NA12878 exome reads (SRR098401) and the GRCh38 reference
# Part 6: download the GIAB HG001 v4.2.1 benchmark (truth VCF + confident BED)
#
# Run from the project folder:   bash scripts/01_download.sh
# Conda environment:             gvp2
# Disk space needed:             ~25 GB
# =============================================================================

set -euo pipefail

# ---- Project folders --------------------------------------------------------
mkdir -p reads reference results truth qc logs

# ---- 1. Run metadata and md5 checksums from ENA ------------------------------
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR098401&result=read_run&fields=sample_alias,instrument_model,library_strategy,library_layout,read_count,base_count,fastq_bytes&format=tsv" > logs/SRR098401_metadata.tsv
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR098401&result=read_run&fields=fastq_md5&format=tsv" > logs/SRR098401_md5_ena.tsv

# ---- 2. Exome reads: paired-end FASTQ (~15 GB) -------------------------------
wget -c -P reads https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098401/SRR098401_1.fastq.gz
wget -c -P reads https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098401/SRR098401_2.fastq.gz

# Integrity check: compare these sums with logs/SRR098401_md5_ena.tsv
md5sum reads/SRR098401_1.fastq.gz reads/SRR098401_2.fastq.gz | tee logs/SRR098401_md5_local.txt

# ---- 3. Reference: GRCh38 no-alt analysis set (NCBI, "chr" names) ------------
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
wget -c -P reference ${NCBI}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget -c -P reference ${NCBI}/md5checksums.txt

# Integrity check: the two sums printed below must be identical
grep "no_alt_analysis_set.fna.gz" reference/md5checksums.txt
md5sum reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Decompress to a short working name; the original .gz is kept
gunzip -c reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > reference/GRCh38_no_alt.fa

# ---- 4. GIAB HG001 benchmark, GRCh38, v4.2.1 (chr1-22) -----------------------
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
wget -c -P truth ${GIAB}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c -P truth ${GIAB}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -c -P truth ${GIAB}/HG001_GRCh38_1_22_v4.2.1_benchmark.bed

echo "Downloads finished."
