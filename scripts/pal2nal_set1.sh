#!/bin/bash

#SBATCH --array=1-773
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/pal2nal_set1_%A_%a.out
#SBATCH --error=log/pal2nal_set1_%A_%a.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/pal2nal.v14:${PATH}

seq=$(cat scripts/busco_ids_set1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

pal2nal.pl analyses/phylogenetics/set1/alignments/single/${seq}_aln.faa \
analyses/phylogenetics/set1/seqs/${seq}.fna -codontable 11 -output fasta > \
analyses/phylogenetics/set1/alignments/single/${seq}_aln.fna
