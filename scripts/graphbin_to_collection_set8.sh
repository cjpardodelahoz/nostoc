#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/graphbin_to_collection_set8_%A_%a.err
#SBATCH --output=log/graphbin_to_collection_set8_%A_%a.out
#SBATCH --partition=scavenger

# Sample ID variable
S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Bin paths variable
S_bins=$(ls -d analyses/genome_qc/set8/delta_contigs/${S}_*)
# Get delta contigs for each cyano bin from a sample
i=1
for bin in ${S_bins} ; do
 cat ${bin}/set8_bandage.csv | sed '1d' | sed "s/$/${i}/" | sed 's/,/\t/' > \
  analyses/bins/anvio/${S}/delta_contigs_${i}
 ((i++))
done
# Merge delta contigs for the sample in anvio collection format
cat analyses/bins/anvio/${S}/delta_contigs_* > \
 analyses/bins/anvio/${S}/graphbin_delta_collection.txt
cat analyses/bins/anvio/${S}/delta_contigs_* | sed 's/\t/_split_00001\t/' | \
 sed '1 i\split\tbin' > \
 analyses/bins/anvio/${S}/graphbin_delta_collection_with_header_splitnames.txt
