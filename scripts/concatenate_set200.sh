#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_set200.out
#SBATCH --error=log/concatenate_set200.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

mkdir -p analyses/phylogenetics/set200/alignments/concat

# Concat alignments with 1 (all genes) and 1648 (partitioned by gene) partitions

# Make variable with paths to alignments
aln_dir=analyses/phylogenetics/set200/alignments/single/
aln_paths=$(cat scripts/busco_ids_set200 | sed 's/$/_ng.faa/' | sed "s|^|${aln_dir}|")
# Concatenate 1648 loci
AMAS.py concat -i ${aln_paths} \
 -f fasta -d aa -p analyses/phylogenetics/set200/alignments/concat/gene_partition_ng_aa \
 --concat-out analyses/phylogenetics/set200/alignments/concat/concat_ng_1_part.phy \
 --out-format phylip
# Split the alignments by gene, which will now include all taxa
AMAS.py split -l analyses/phylogenetics/set200/alignments/concat/gene_partition_ng_aa \
 --out-format phylip \
 -i analyses/phylogenetics/set200/alignments/concat/concat_ng_1_part.phy \
 -d aa -f phylip
# Move the phylip splitted alignments to their own directory
mkdir -p analyses/phylogenetics/set200/alignments/concat/splitted
mv analyses/phylogenetics/set200/alignments/concat/*out.phy \
 analyses/phylogenetics/set200/alignments/concat/splitted
# Compile all phylip gene alignments into a single file, which will serve as the
# concatenated partitioned file for MCMCTree
cat analyses/phylogenetics/set200/alignments/concat/splitted/*.phy > \
 analyses/phylogenetics/set200/alignments/concat/concat_ng_1648_part.phy
 
# Concat alignment with 4 partitions each with 412 loci

# Split the busco ids file into 4 random sets of 412
ls -p analyses/phylogenetics/set200/alignments/concat/splitted/*phy| sort -R | \
 split --lines 412 --suffix-length 1 \
 - analyses/phylogenetics/set200/alignments/concat/busco_ids_
# Concatenate the genes for each set
for id in {a..d} ; do
 AMAS.py concat -i $(cat analyses/phylogenetics/set200/alignments/concat/busco_ids_${id}) \
 -f phylip -d aa --out-format phylip \
 --concat-out analyses/phylogenetics/set200/alignments/concat/concat_ng_${id}.faa
done
# Concatenate the for partitions in phylip format
cat analyses/phylogenetics/set200/alignments/concat/concat_ng_{a..d}.faa > \
 analyses/phylogenetics/set200/alignments/concat/concat_ng_4_part.phy
 