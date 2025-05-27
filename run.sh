#!/bin/bash
#SBATCH --job-name=rnaseq_nf
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

mkdir -p logs results
module load nextflow
module load singularity

nextflow run main.nf -profile slurm --reads 'raw_fastq/*_{1,2}.fastq.gz' -c nextflow.config
