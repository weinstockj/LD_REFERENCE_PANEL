# LD Reference Panel Generation

This Snakemake pipeline generates a directory structure suitable for use as an LD reference panel for PRSFNN (Polygenic Risk Score Neural Network). The pipeline processes 1000 Genomes Project VCF files to create population-specific genotype files in PLINK format, organized by LD blocks.

## Overview

The pipeline performs the following steps:
1. Converts VCF files to BCF format with quality filtering
2. Extracts population-specific samples based on 1000 Genomes metadata
3. Filters variants by LD blocks (genomic regions in linkage disequilibrium)
4. Applies population-specific MAF filtering (>0.5% by default)
5. Converts filtered data to PLINK format (.bed/.bim/.fam files)
6. Generates variant lists for downstream analysis

## Input Data

- **VCF files**: 1000 Genomes Project variant call files (one per chromosome)
- **LD blocks file**: BED format file defining LD block boundaries (currently supports EUR blocks)
- **Population file**: 1000 Genomes sample metadata with population assignments

## Output Structure

```
output/bcf/
├── {chr}.recalibrated_variants.bcf     # Chromosome-wide BCF files
├── {pop}.samples                       # Population sample lists
└── {block}/                           # LD block directories
    ├── filtered.bcf                   # All samples filtered by block
    ├── filtered_{pop}.bcf             # Population-specific filtered data
    ├── filtered.bed/.bim/.fam         # PLINK format files (all samples)
    ├── filtered_{pop}.bed/.bim/.fam   # PLINK format files (population-specific)
    ├── variant_list.tsv               # Variant site list (all samples)
    └── variant_list_{pop}.tsv         # Variant site list (population-specific)
```

## Configuration

The pipeline reads configuration parameters from `config.yaml`. Key parameters include:
- Input/output file paths
- Population groups to process
- Quality control thresholds

## Usage

1. Configure parameters in `config.yaml`
2. Run the pipeline:
   ```bash
   snakemake --profile slurm --jobs 100
   ```
   or use the provided cluster script:
   ```bash
   bash run_cluster.sh
   ```

## Dependencies

**Required software** (note: other versions may work but are untested):
- PLINK v2.00a5.9LM AVX2 Intel (12 Dec 2023)
- bcftools version 1.15
- Snakemake 8.0.1
- Python packages: pandas, snakemake

**System requirements**:
- HPC cluster with SLURM scheduler (for parallel processing)
- 1000 Genomes ~30x WGS

## Quality Control

The pipeline applies several QC filters:
- Biallelic SNPs only
- PASS filter from VCF
- Minor allele count ≥ 5 (whole dataset)
- Missing data rate < 15%
- Population-specific MAF > 0.5% (configurable)

## Notes

- Currently optimized for European (EUR) LD blocks
- No X chromosome processing (autosomal chromosomes only)  
- Designed for use with SLURM job scheduler
- Output files are indexed for efficient random access
