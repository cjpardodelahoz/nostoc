#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/get_16s_from_ffn_set103_%A_%a.err
#SBATCH --output=log/get_16s_from_ffn_set103_%A_%a.out
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/:$PATH

bin=$(cat scripts/genome_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)

# .ffn file output from Prokka annotation
ffn_file=analyses/cyano_genomes_annotation/set103/${bin}/${bin}.ffn
# The full header of the first 16 sequence in the annotation file
seq_id=$(grep "16S ribosomal RNA" ${ffn_file} | head -1)
# If the 16S is present..
if [[ -n ${seq_id} ]]; then
 # Get the seq ID as interpreted by seqkit for the 16 sequence
 seq_pattern=$(echo ${seq_id} | head -1 | sed 's/>//' | sed 's/ 16S ribosomal RNA//')
 # Path to output file
 outfile=analyses/cyano_genomes_annotation/set103/${bin}/16s.fas
 # Get sequence and rename it with bin name
 seqkit grep -p ${seq_pattern} ${ffn_file} --immediate-output -r | \
 sed 's/ 16S ribosomal RNA//' | sed "s/${seq_pattern}/${bin}/" > ${outfile}
# If the 16s is not present, print the name of the bin
else
 echo ${bin}
fi
