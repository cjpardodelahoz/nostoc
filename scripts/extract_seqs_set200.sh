#!/bin/bash

#SBATCH --array=1-1648
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/extract_seqs_set200_%A_%a.err
#SBATCH --output=log/extract_seqs_set200_%A_%a.out
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/:$PATH

busco=$(cat scripts/busco_ids_set200 | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Extract the set200 taxa from the original Pardo De la Hoz sequence files
seqkit grep -f scripts/genome_ids_set200 \
 analyses/phylogenetics/set200/pardodelahoz2023_seqs/${busco}*.faa \
 -o analyses/phylogenetics/set200/seqs/${busco}.faa