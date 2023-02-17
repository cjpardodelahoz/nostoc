#!/bin/bash

#SBATCH --array=1-58
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/fastqc_6240_6526_%A_%a.out
#SBATCH --error=log/fastqc_6240_6526_%A_%a.err
#SBATCH --partition=common

module load FastQC/0.11.7

S=$(cat scripts/sample_ids_6240_6526 | sed -n ${SLURM_ARRAY_TASK_ID}p)

fastqc --outdir analyses/reads/${S}/ -t 4 analyses/reads/${S}/*all.fastq.gz
