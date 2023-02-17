#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/spades_vecherskii_%A_%a.out
#SBATCH --error=log/spades_vecherskii_%A_%a.err
#SBATCH --partition=common

module load Python/3.8.1
module load SPAdes/3.14.1

S=$(cat scripts/sample_ids_vecherskii | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p /work/cjp47/nostoc/assemblies/${S}
ln -s /work/cjp47/nostoc/assemblies/${S} analyses/assemblies/${S}

spades.py --meta --pe1-1 analyses/reads/${S}/${S}_R1_paired.fq.gz \
--pe1-2 analyses/reads/${S}/${S}_R2_paired.fq.gz -k 25,35,55,75,81 -t 16 \
-o analyses/assemblies/${S}
