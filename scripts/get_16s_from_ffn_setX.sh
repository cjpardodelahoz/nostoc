#!/bin/bash

#SBATCH --array=1
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/get_16s_from_ffn_%A_%a.err
#SBATCH --output=log/get_16s_from_ffn_%A_%a.out
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/:$PATH

bin=$(cat scripts/genome_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

# .ffn file output from Prokka annotation
ffn_file=analyses/cyano_genomes/set8/${bin}_annotation/${bin}.ffn
# The full header of the 16 sequence in the annotation file
seq_id=$(grep "16S ribosomal RNA" ${ffn_file})
# The seq ID as interpreted by seqkit for the 16 sequence
seq_pattern=$(echo ${seq_id} | sed 's/>//' | sed 's/ 16S ribosomal RNA//')
# Path to output file
outfile=analyses/cyano_genomes/set8/${bin}_annotation/16s.fas

seqkit grep -p ${seq_pattern} ${ffn_file} --immediate-output -r | \
seqkit replace -p ${seq_id} -r ${bin} > ${outfile}
