#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH --error=log/contigs_to_graph_set7_%A_%a.err
#SBATCH --output=log/contigs_to_graph_set7_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load Python/2.7.11
export PATH=/hpc/group/bio1/carlos/apps/SPAdes-Contig-Graph:$PATH

spades_contig_graph.py -l analyses/assemblies/${S}/assembly_graph.fastg \
analyses/assemblies/${S}/contigs.fasta analyses/assemblies/${S}/contigs.paths \
analyses/assemblies/${S}/contigs.fastg
