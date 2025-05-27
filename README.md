
This is a modular and reproducible RNA-Seq pipeline written in [Nextflow](https://www.nextflow.io/). It performs quality control, alignment, gene-level quantification, and transcript-level quantification for paired-end RNA-Seq FASTQ files.

## Downloading raw FASTQ files

First, make sure you have `parallel-fastq-dump` installed (e.g. via conda):

```bash
conda install -c bioconda parallel-fastq-dump
```
Then for each SRR ID, 
```bash
parallel-fastq-dump \
  --sra-id SRR1234567     # SRA accession
  --threads 16            # number of CPU threads to use
  --split-files           # split into *_1.fastq and *_2.fastq
  --gzip                  # compress output with gzip
  --outdir raw_fastq/     # output directory
```

## Requirements
Nextflow (tested on v24.10.5)
Singularity (tested onv3.11.4)
HPC with SLURM

####

# RNA-Seq Pipeline using Nextflow

A containerized RNA-seq analysis pipeline built with Nextflow and Singularity, featuring dual quantification using STAR alignment and Salmon pseudo-alignment.

## Overview

This pipeline processes paired-end RNA-seq data through quality control, alignment, and quantification steps using tools packaged in Singularity containers for reproducible execution across computing environments.

### Features

- **Quality Control**: FastQC analysis of raw reads
- **Dual Quantification**: 
  - Genome alignment with STAR followed by featureCounts
  - Transcript-level quantification with Salmon
- **Containerized**: All tools run in Singularity containers

## Pipeline Steps

1. **FastQC** - Quality assessment of raw FASTQ files
2. **STAR Index** - Build genome index for alignment
3. **Salmon Index** - Build transcriptome index for quantification
4. **STAR Alignment** - Map reads to reference genome
5. **featureCounts** - Gene-level read counting from alignments
6. **Salmon Quantification** - Transcript-level abundance estimation

## Requirements

### Software Dependencies
- Nextflow (tested on v24.10.5)
- Singularity (tested onv3.11.4)
- HPC with SLURM

### Reference Files
- Reference genome FASTA
- Transcriptome FASTA
- GTF annotation file

## Configuration

Edit `nextflow.config` to specify your reference files:

```groovy
params {
    reads = "raw_fastq/*_{1,2}.fastq.gz"
    outdir = "results"
    gtf = "/path/to/gencode.v47.annotation.gtf"
    genome_fasta = "/path/to/GRCh38.p14.genome.fa"
    transcriptome_fasta = "/path/to/gencode.v47.transcripts.fa"
}
```

## Usage

### Local Execution
```bash
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz'
```

### SLURM Execution
```bash
# Submit via SLURM job script
sbatch run.sh

# Or run directly with SLURM profile
nextflow run main.nf -profile slurm --reads 'data/*_{1,2}.fastq.gz'
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Path to paired-end FASTQ files | `"path/to/fastqs/*_{1,2}.fastq.gz"` |
| `--outdir` | Output directory | `"results"` |
| `--gtf` | GTF annotation file | Required |
| `--genome_fasta` | Reference genome FASTA | Required |
| `--transcriptome_fasta` | Transcriptome FASTA | Required |

## Output Structure

```
results/
├── fastqc/           # FastQC reports
├── align/            # STAR alignment BAM files
├── star_index/       # STAR genome index
├── salmon_index/     # Salmon transcriptome index
└── salmon_quant/     # Salmon quantification results
```

## Container Images

All tools are executed within Singularity containers pulled from Biocontainers:

- **FastQC**: `quay.io/biocontainers/fastqc:0.11.9--0`
- **STAR**: `quay.io/biocontainers/star:2.7.11b--h5ca1c30_6`
- **Subread**: `quay.io/biocontainers/subread:2.1.1--h577a1d6_0`
- **Salmon**: `quay.io/biocontainers/salmon:1.10.3--h45fbf2d_4`

Resource requirements can be adjusted in the `nextflow.config` file or via command-line parameters.

## Citations

- Nextflow: Di Tommaso et al. (2017) Nature Biotechnology
- STAR: Dobin et al. (2013) Bioinformatics  
- Salmon: Patro et al. (2017) Nature Methods
- FastQC: Andrews (2010) Babraham Bioinformatics
