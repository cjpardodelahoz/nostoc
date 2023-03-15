#!/bin/bash

#SBATCH --array=1-124
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ncbi_fcs_set103p_%A_%a.out
#SBATCH --error=log/ncbi_fcs_set103p_%A_%a.err
#SBATCH --partition=scavenger

# FCS script on path and sif
export PATH=/hpc/group/bio1/carlos/apps/fcsadaptor:${PATH}
fcs_sif="/hpc/group/bio1/carlos/apps/fcsadaptor/fcs-adaptor.sif"
# Variable with chromosome files
genome=$(cat misc_files/genome_ids_set10 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Make output directory
mkdir -p analyses/genome_qc/set103p/fcs/${genome%.fa}_plasmid.fa 
# Run FCS
run_fcsadaptor.sh --fasta-input analyses/cyano_genomes/set103p/${genome%.fa}_plasmid.fa  \
 --output-dir analyses/genome_qc/set103p/fcs/${genome%.fa}_plasmid.fa \
 --prok --container-engine singularity \
 --image ${fcs_sif}

