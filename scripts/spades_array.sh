#!/bin/bash

#SBATCH --array=1-32
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 14 # number of threads per process
#SBATCH --output=log/spades_%A_%a.out
#SBATCH --error=log/spades_%A_%a.err
#SBATCH --partition=common

module load Python/3.8.1
module load SPAdes/3.14.1

S=$(cat scripts/sample_ids | sed -n ${SLURM_ARRAY_TASK_ID}p)

spades.py --meta --pe1-1 analyses/reads/${S}/${S}_R1_paired.fq.gz \
--pe1-2 analyses/reads/${S}/${S}_R2_paired.fq.gz -k 55,75,95 -t 14 \
-m 112 -o analyses/assemblies/${S}
