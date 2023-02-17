#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=4G
#SBATCH -c 8
#SBATCH --error=log/graphbin_set7_%A_%a.err
#SBATCH --output=log/graphbin_set7_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)

source $(conda info --base)/etc/profile.d/conda.sh
conda activate graphbin2
export PATH=/hpc/group/bio1/carlos/apps/GraphBin2:${PATH}

graphbin2 --nthreads 8 --assembler spades --contigs analyses/assemblies/${S}/contigs.fasta \
--graph analyses/assemblies/${S}/assembly_graph_after_simplification.gfa \
--paths analyses/assemblies/${S}/contigs.paths --binned analyses/bins/graphbin/binned_csv/${S}_contig_labels.csv \
--output analyses/bins/graphbin/${S}