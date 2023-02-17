#!/bin/bash

#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=graph.err
#SBATCH --output=graph.out
#SBATCH --partition=scavenger

module load Python/2.7.11
export PATH=/hpc/group/bio1/carlos/apps/SPAdes-Contig-Graph:$PATH

spades_contig_graph.py -l analyses/assemblies/P6465/assembly_graph.fastg analyses/assemblies/P6465/contigs.fasta analyses/assemblies/P6465/contigs.paths contigs.fastg
