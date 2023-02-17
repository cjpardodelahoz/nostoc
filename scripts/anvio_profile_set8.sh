#!/bin/bash

#SBATCH --array=22
#SBATCH --mem-per-cpu=38G  # adjust as needed
#SBATCH -c 18 # number of threads per process
#SBATCH --output=log/anvio_profile_set8_%A_%a.out
#SBATCH --error=log/anvio_profile_set8_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

rm -rf analyses/bins/anvio/${S}/profile_1

anvi-profile -c analyses/bins/anvio/${S}/contigs.db \
 -i analyses/assemblies/${S}/${S}_sorted.bam \
 -o analyses/bins/anvio/${S}/profile_1 \
 --sample-name s_${S} \
 --cluster-contigs \
 --min-contig-length 90 \
 --num-threads 18
 
 
 
 