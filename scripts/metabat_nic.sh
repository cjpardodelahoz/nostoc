#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=25G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/metabat_nic_%A_%a.out
#SBATCH --error=log/metabat_nic_%A_%a.err
#SBATCH --partition=scavenger

module load MetaBAT/2020-10-09

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

/opt/apps/rhel7/metabat-2020-10-09/bin/metabat2 -t 8 -i analyses/assemblies/${S}/contigs.fasta \
-a analyses/assemblies/${S}/${S}_assembly_depths.txt -o analyses/bins/metabat/${S}/${S}_bin
