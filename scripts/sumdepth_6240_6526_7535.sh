#!/bin/bash

#SBATCH --array=90
#SBATCH --mem-per-cpu=5G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/sumdepth_6240_6526_7535_%A_%a.out
#SBATCH --error=log/sumdepth_6240_6526_7535_%A_%a.err
#SBATCH --partition=scavenger

module load MetaBAT/2020-10-09

S=$(cat scripts/sample_ids_6240_6526_7535 | sed -n ${SLURM_ARRAY_TASK_ID}p)

/opt/apps/rhel7/metabat-2020-10-09/bin/jgi_summarize_bam_contig_depths \
--outputDepth analyses/assemblies/${S}/${S}_assembly_depths.txt \
analyses/assemblies/${S}/${S}_sorted.bam
