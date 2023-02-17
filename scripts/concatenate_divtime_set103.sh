#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_divtime_set103.out
#SBATCH --error=log/concatenate_divtime_set103.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

mkdir -p analyses/phylogenetics/set200/alignments/concat

# Concat alignments with 1 (all 1519 genes) partitions

# Make variable with paths to filtered alignments
aln_dir=analyses/phylogenetics/set103/alignments/single/
aln_paths=$(cat scripts/busco_ids_filtered_set103 | sed 's/$/_ng.fna/' |  sed '1i 16s_aln.fas' | sed '1i trnl_aln.fas' | sed "s|^|${aln_dir}|")
# Concatenate 1519 loci
AMAS.py concat -i ${aln_paths} \
 -f fasta -d dna -p analyses/phylogenetics/set103/alignments/concat/gene_partition_ng_na.phy \
 --concat-out analyses/phylogenetics/set103/alignments/concat/concat_ng_1_part.phy \
 --out-format phylip
# Split the alignments by gene, which will now include all taxa
AMAS.py split -l analyses/phylogenetics/set103/alignments/concat/gene_partition_ng_na.phy \
 --out-format phylip \
 -i analyses/phylogenetics/set103/alignments/concat/concat_ng_1_part.phy \
 -d dna -f phylip
# Move the phylip splitted alignments to their own directory
mkdir -p analyses/phylogenetics/set103/alignments/concat/splitted
mv analyses/phylogenetics/set103/alignments/concat/*out.phy \
 analyses/phylogenetics/set103/alignments/concat/splitted
 
# Concat alignment with 4 partitions each with 217 loci

# Split the busco ids file into 4 random sets of 380
ls -p analyses/phylogenetics/set103/alignments/concat/splitted/*phy| sort -R | \
 split --lines 380 --suffix-length 1 \
 - analyses/phylogenetics/set103/alignments/concat/busco_ids_
# Concatenate the genes for each set
for id in {a..d} ; do
 AMAS.py concat -i $(cat analyses/phylogenetics/set103/alignments/concat/busco_ids_${id}) \
 -f phylip -d dna --out-format phylip \
 --concat-out analyses/phylogenetics/set103/alignments/concat/concat_ng_${id}.fna
done
# Concatenate the four partitions in phylip format
cat analyses/phylogenetics/set103/alignments/concat/concat_ng_{a..d}.fna > \
 analyses/phylogenetics/set103/alignments/concat/concat_ng_4_part.phy
 