# Genomics Formats Quickstart

A minimal Python project to inspect common genomics file formats: **FASTQ**, **BAM**, and **VCF**.

The goal is to get hands-on with the exact data formats used in genomics and pathology pipelines (NGS, variant calling, spatial transcriptomics) without needing a full bioinformatics workflow.

## Data

All sample files are small, real public datasets from the 1000 Genomes Project:

- `NA12878_R1_subset.fastq` + `NA12878_R2_subset.fastq` — 100 paired-end Illumina reads from sample NA12878 ( chromosome 20 region, ERR001268).
- `sample.bam` / `sample.bam.bai` — small indexed BAM alignment subset.
- `basic.vcf` — single-sample VCF from 1000 Genomes.
- `basic_multisample.vcf` — multi-sample VCF (1,233 samples).
- `indexed_tbi.vcf.gz` / `indexed_tbi.vcf.gz.tbi` — bgzip-compressed, tabix-indexed VCF.

## Requirements

```bash
pip install -r requirements.txt
```

- Python 3.10+
- `biopython` — FASTQ parsing
- `pysam` — BAM/SAM parsing
- `cyvcf2` — VCF/BCF parsing

## Run

```bash
python scripts/inspect_genomics_formats.py
```

Outputs are written to `outputs/`:
- JSON summaries for each format
- CSV tables for distributions (read lengths, mapping quality, variants per chromosome, filters, etc.)

## What the script demonstrates

- **FASTQ**: read counts, read length distribution, mean Phred quality per file.
- **BAM**: mapped/unmapped counts, reference distribution, mapping quality distribution, CIGAR operation counts, insert-size distribution.
- **VCF**: variant counts, SNP/indel split, filter counts, chromosome distribution, allele counts, genotype counts (for multi-sample files).

## Next steps

- Point the script at larger FASTQ/BAM/VCF files from TCGA, 1000 Genomes, or HEST-1k.
- Connect a WSI TIFF workflow with `tifffile` / `openslide-python` to bridge into pathology.
- Run `nf-core/sarek` or a GATK/DeepVariant quick-start to see the full FASTQ → BAM → VCF pipeline end-to-end.

## License

MIT
