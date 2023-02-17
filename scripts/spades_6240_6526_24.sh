#!/bin/bash

#SBATCH --array=24
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/spades_6240_6526_%A_%a.out
#SBATCH --error=log/spades_6240_6526_%A_%a.err
#SBATCH --partition=common

module load Python/3.8.1
module load SPAdes/3.14.1

S=$(cat scripts/sample_ids_6240_6526 | sed -n ${SLURM_ARRAY_TASK_ID}p)

spades.py -o analyses/assemblies/${S} --continue
