#!/bin/bash

#SBATCH --array=1-1898
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/trimal_ng_set103_%A_%a.out
#SBATCH --error=log/trimal_ng_set_set103_%A_%a.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/trimAl/source:${PATH}

seq=$(cat scripts/busco_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)

trimal -in analyses/phylogenetics/set103/alignments/single/${seq}_aln.faa \
-out analyses/phylogenetics/set103/alignments/single/${seq}_ng.faa \
-fasta -nogaps

trimal -in analyses/phylogenetics/set103/alignments/single/${seq}_aln.fna \
-out analyses/phylogenetics/set103/alignments/single/${seq}_ng.fna \
-fasta -nogaps
