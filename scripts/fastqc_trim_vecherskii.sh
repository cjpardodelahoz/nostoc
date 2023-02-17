#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/fastqc_trim_vecherskii_%A_%a.out
#SBATCH --error=log/fastqc_trim_vecherskii_%A_%a.err
#SBATCH --partition=scavenger

module load FastQC/0.11.7

sample=$(cat scripts/sample_ids_vecherskii | sed -n ${SLURM_ARRAY_TASK_ID}p)

fastqc --outdir analyses/reads/${sample}/ -t 4 analyses/reads/${sample}/*_paired.fq.gz
