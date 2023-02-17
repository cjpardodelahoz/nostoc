#!/bin/bash

#SBATCH --array=1-129
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/get_delta_contigs_set8_%A_%a.out
#SBATCH --error=log/get_delta_contigs_set8_%A_%a.err
#SBATCH --partition=scavenger

bin=$(cat scripts/genome_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p analyses/genome_qc/set8/delta_contigs/${bin}
# Get list of contigs for a given genome in set7 and set8
grep ">" analyses/cyano_genomes/set7/${bin} | sed 's/>//' > \
analyses/genome_qc/set8/delta_contigs/${bin}/set7_all_contigs
grep ">" analyses/cyano_genomes/set8/${bin} | sed 's/>//' > \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_all_contigs
# Get added contigs csv
grep -vxF -f analyses/genome_qc/set8/delta_contigs/${bin}/set7_all_contigs \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_all_contigs | \
sed 's/$/,green/' > \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_added_contigs
# Get removed contigs
grep -vxF -f analyses/genome_qc/set8/delta_contigs/${bin}/set8_all_contigs \
analyses/genome_qc/set8/delta_contigs/${bin}/set7_all_contigs |
sed 's/$/,red/' > \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_removed_contigs
# Merge added and removed contigs
cat analyses/genome_qc/set8/delta_contigs/${bin}/set8_added_contigs \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_removed_contigs > \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_delta
# Get unchaged contigs
cat analyses/genome_qc/set8/delta_contigs/${bin}/set8_delta | sed 's/,.*//' > \
 analyses/genome_qc/set8/delta_contigs/${bin}/set8_delta_labels
grep -vxF -f analyses/genome_qc/set8/delta_contigs/${bin}/set8_delta_labels \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_all_contigs | \
sed 's/$/,black/' > \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_unchanged_contigs
# Merge delta and unchanged contigs
cat analyses/genome_qc/set8/delta_contigs/${bin}/set8_delta \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_unchanged_contigs | \
sed '1 i\name,color'> \
analyses/genome_qc/set8/delta_contigs/${bin}/set8_bandage.csv




