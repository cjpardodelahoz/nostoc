#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/get_codon_partition_concat_set103.out
#SBATCH --error=log/get_codon_partition_concat_set103.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

partition_dir=analyses/phylogenetics/set103/alignments/concat
# generate codon partition file from gene partition file - ng
for line in $(sed 's/\ /./g' ${partition_dir}/gene_partition_ng_na) ; do
# get the locus name from the partition definition line
	locus=$(echo ${line} | grep -e "DNA" | awk -F'[.]' '{print $2}')
# get first and last position of the gene in the concatenated alignment
	start=$(echo ${line} | grep -e "DNA" | awk -F'[.]' '{print $4}' | awk -F'[-]' '{print $1}')
	fin=$(echo ${line} | grep -e "DNA" | awk -F'[.]' '{print $4}' | awk -F'[-]' '{print $2}')
# write line for 1st, 2nd and 3rd codon position
	printf "DNA, 1st_${locus} = ${start}-${fin}\n"
	printf "DNA, 2nd_${locus} = $(($start + 1))-${fin}\n"
	printf "DNA, 3rd_${locus} = $(($start + 2))-${fin}\n"
done | sed -e 's/$/\\3/' > ${partition_dir}/codon_partition_concat_ng_na