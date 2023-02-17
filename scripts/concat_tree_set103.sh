#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -n 8 # number of processes
#SBATCH -c 12 # number of threads per process
#SBATCH --output=log/concat_tree_set103.out
#SBATCH --error=log/concat_tree_set103.err
#SBATCH --partition=common

module load IQ-TREE/1.6.12-MPI

mkdir -p analyses/phylogenetics/set103/trees/concat

mpirun -np 8 iqtree-mpi -nt 12 -s analyses/phylogenetics/set103/alignments/concat/concat_ng.fna \
 -spp analyses/phylogenetics/set103/trees/concat/concat_pf_ng_na.best_scheme.nex \
 -pre analyses/phylogenetics/set103/trees/concat/concat_tree_ng_na \
 -bb 1000 -safe 
