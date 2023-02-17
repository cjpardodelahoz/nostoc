#!/bin/bash

#SBATCH --array=1-95
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_genomes2020_resequencing2021_%A_%a.err
#SBATCH --output=log/merge_reads_genomes2020_resequencing2021_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_genomes2020_resequencing2021 | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Copy reads into analyses directory
# If one of the sample files exists in the directory, copy all of them
if [ -e data/reads/nic_reads/resequencing2021/6126-${S}_*L001_R1_001.fastq.gz ]
then
 cp data/reads/nic_reads/resequencing2021/6126-${S}_*.gz analyses/reads/
fi

if [ -e data/reads/nic_reads/genomes2020/6126-${S}_*L001_R1_001.fastq.gz ]
then
 cp data/reads/nic_reads/genomes2020/6126-${S}_*.gz analyses/reads/
fi
# Make sample specific directory
mkdir analyses/reads/${S}
# Merge reads and sort into folders
cat analyses/reads/6126-${S}_*R1* > analyses/reads/${S}/${S}_R1_all.fastq.gz
cat analyses/reads/6126-${S}_*R2* > analyses/reads/${S}/${S}_R2_all.fastq.gz
# Delete copied read files
rm analyses/reads/6126-${S}_*R1*
rm analyses/reads/6126-${S}_*R2*
