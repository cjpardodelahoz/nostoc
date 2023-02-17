#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/get_codon_partition_set103.out
#SBATCH --error=log/get_codon_partition_set103.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

aln_path=analyses/phylogenetics/set103/alignments/single
# initiate partition file using AMAS hack
for locus in $(cd ${aln_path} && ls *ng.fna | awk -F'[.]' '{print $1}') ; do
	AMAS.py concat -i ${aln_path}/${locus}.fna -f fasta -d dna -p ${aln_path}/${locus}_part --part-format raxml
# modify AMAS file to get codon parition file
	grep -e "DNA" ${aln_path}/${locus}_part | sed 's/p1/p2/' | sed 's/1-/2-/' >> ${aln_path}/${locus}_part 
	grep -e "DNA" ${aln_path}/${locus}_part | sed -n 1p | sed 's/p1/p3/' | sed 's/1-/3-/' >> ${aln_path}/${locus}_part
	cat ${aln_path}/${locus}_part | sed -e 's/$/\\3/' > ${aln_path}/${locus}_codon_partition
done
# remove intermediate partition files
rm ${aln_path}/*_part
rm concatenated.out