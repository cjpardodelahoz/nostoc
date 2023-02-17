#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ln_assemblies_nic_%A_%a.out
#SBATCH --error=log/ln_assemblies_nic_%A_%a.err
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

#cp -r analyses/assemblies/${S} /work/cjp47/nostoc/assemblies
rm -rf analyses/assemblies/${S}
ln -s /work/cjp47/nostoc/assemblies/${S} analyses/assemblies/${S}
