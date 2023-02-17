#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -n 4 # number of processes
#SBATCH -c 12 # number of threads per process
#SBATCH --output=log/concat_ng_set1_na1.out
#SBATCH --error=log/concat_ng_set1_na1.err
#SBATCH --partition=common

module load IQ-TREE/1.6.12-MPI

mpirun -np 4 iqtree-mpi -nt 12 -s analyses/phylogenetics/set1/alignments/concat/concat_ng.fna \
-m GTR+G -bb 1000 -pre analyses/phylogenetics/set1/trees/concat/concat_ng_na1
