# WES Germline Variant Analysis of NA12878

End-to-end whole exome sequencing (WES) germline variant analysis of the reference sample NA12878 (GIAB HG001): from raw reads to annotated, prioritized variants, benchmarked against the Genome in a Bottle (GIAB) v4.2.1 truth set with hap.py.

**Headline result** (PASS calls, final filter `QUAL<30 || DP<10 || MQ<40`, 80.1 Mb evaluation region):

| Variant type | Precision | Recall | F1 |
|---|---|---|---|
| SNPs | 0.986 | 0.947 | **0.966** |
| Indels | 0.859 | 0.688 | **0.764** |

---

## Data

| Item | Details |
|---|---|
| Sample | NA12878 / HG001 (CEPH/Utah, B-lymphoblastoid cell line GM12878) |
| Reads | ENA/SRA run [SRR098401](https://www.ebi.ac.uk/ena/browser/view/SRR098401): Illumina HiSeq 2000, 2 x 76 bp paired-end, 114,059,565 read pairs |
| Capture | Agilent SureSelect v2 exome |
| Reference | GRCh38 no-alt analysis set (NCBI, "chr" names) |
| Truth set | GIAB HG001 GRCh38 v4.2.1 benchmark VCF + high-confidence BED (chr1-22) |

Raw data are not stored in this repository; `scripts/01_download.sh` downloads everything and checks md5 checksums.

## Workflow

| Step | Script | Tools | Main output |
|---|---|---|---|
| 1. Download | `01_download.sh` | wget, curl, md5sum | FASTQ, reference, GIAB files |
| 2. QC and alignment | `02_qc_align.sh` | FastQC, bwa mem, samtools | duplicate-marked BAM |
| 3. Calling and filtering | `03_call_filter.sh` | bcftools mpileup/call/filter/norm | normalized VCF, PASS VCF |
| 4. Annotation and prioritization | `04_annotate.sh` | Ensembl VEP 116, bcftools +split-vep | three tier tables |
| 5. Evaluation region | `05_eval_regions.sh` | mosdepth, bedtools, bcftools | callable high-confidence BED |
| 6. Benchmark | `06_benchmark.sh` | hap.py | precision / recall / F1 |

## Results

### Alignment
- 93.66% of reads mapped, 91.25% properly paired.
- 17,958,217 duplicate reads (8.4%) marked with `samtools markdup`.

### Variant calling
- 2,926,174 raw calls; 101,713 PASS after soft filtering (most LowQual calls are low-depth off-target reads).
- PASS Ts/Tv = 2.33.

### Annotation and prioritization (VEP 116, one consequence per variant)
- **Tier 1** (ClinVar pathogenic/likely pathogenic, MAX_AF < 1%): 3 variants
- **Tier 2** (HIGH impact, MAX_AF < 1%, MANE Select): 12 variants
- **Tier 3** (SIFT deleterious + PolyPhen probably damaging, MAX_AF < 1%, VAF >= 0.25): 31 variants

The three Tier 1 variants (TPP1 c.509-1G>A, PAH p.Tyr414Cys, PLG p.Lys38Glu) are all heterozygous (VAF 0.46-0.54) in genes for autosomal recessive conditions, consistent with carrier status in a healthy individual. Placing exon `rank` before `appris` in VEP's `--pick_order` was needed to report TPP1 c.509-1G>A as a splice-acceptor variant rather than an upstream variant of a neighbouring gene.

### Benchmark against GIAB (hap.py)
The evaluation region is the overlap of **callable** regions (>= 10x depth, MAPQ >= 20, from the BAM with mosdepth: 89.0 Mb) and GIAB **high-confidence** regions: **80.1 Mb**, containing 81,891 truth variants.

| Filter (PASS) | SNP precision | SNP recall | SNP F1 | Indel precision | Indel recall | Indel F1 |
|---|---|---|---|---|---|---|
| MQ<50 (original) | 0.988 | 0.924 | 0.955 | 0.859 | 0.682 | 0.760 |
| **MQ<40 (tuned)** | **0.986** | **0.947** | **0.966** | **0.859** | **0.688** | **0.764** |

**Filter tuning.** The benchmark showed the original `MQ<50` threshold was too strict: it removed 5,043 true variants in the evaluation region, 63% of them through the MQ condition. Following the GATK hard-filter recommendation (`MQ < 40`), one re-test rescued 1,718 true SNPs at the cost of 191 new false positives (about 9:1), raising SNP F1 from 0.955 to 0.966. Further tuning against the same truth set was avoided to prevent overfitting to the benchmark.

**Error analysis.** Inspecting individual errors showed four main types:
- **Sampling misses** at low depth (e.g. a true heterozygote seen in only 2 of 8 reads).
- **Filter losses**: correct calls labelled LowQual.
- **Indel genotype errors**: right indel, wrong genotype (half of the indel false positives).
- **Paralog artifacts**: clusters of false heterozygous SNPs in duplicated genes (about 40x false-positive enrichment at chr1:1.63-1.73 Mb), where the alternate-allele reads have low mapping quality (bcftools `MQBZ` = -6.5).

## Software

| Environment | Key tools | File |
|---|---|---|
| `gvp2` | bwa 0.7.19, samtools 1.24, bcftools 1.24, FastQC 0.13.0, mosdepth 0.3.14, bedtools 2.31.1 | `envs/gvp2.yml` |
| `gvp2-vep` | Ensembl VEP 116.3, bcftools 1.21 | `envs/gvp2-vep.yml` |
| `gvp2-happy` | hap.py 0.3.15 (Python 2.7) | `envs/gvp2-happy.yml` |

## How to reproduce

```bash
git clone https://github.com/Hakan-Cam/wes-germline-variant-analysis-na12878.git
cd wes-germline-variant-analysis-na12878

conda env create -f envs/gvp2.yml
conda env create -f envs/gvp2-vep.yml
conda env create -f envs/gvp2-happy.yml

conda activate gvp2
bash scripts/01_download.sh
bash scripts/02_qc_align.sh
bash scripts/03_call_filter.sh

conda activate gvp2-vep
bash scripts/04_annotate.sh

conda activate gvp2
bash scripts/05_eval_regions.sh

conda activate gvp2-happy
bash scripts/06_benchmark.sh
```

Requirements: Linux or WSL2, about 100 GB free disk space, 16+ GB RAM. Alignment takes several hours.

## Repository structure

```
envs/       conda environment files
scripts/    pipeline scripts, run in order 01-06
logs/       QC, alignment, calling, annotation and benchmark logs
qc/         FastQC reports
results/    tier tables (Part 5) and hap.py summaries (results/benchmark/)
```

## Limitations and future work
- Indel recall is limited by 76 bp reads and bcftools' conservative indel calling; a modern caller (e.g. DeepVariant or GATK HaplotypeCaller) would be a natural comparison.
- A mapping-quality bias filter (bcftools `MQBZ`, analogous to GATK `MQRankSum`) could remove paralog artifacts.
- Benchmarking covers chr1-22 only, as defined by the GIAB v4.2.1 truth set.

## Author
Hakan Cam, Ph.D. ([GitHub: Hakan-Cam](https://github.com/Hakan-Cam))
