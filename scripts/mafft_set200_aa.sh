#!/bin/bash

#SBATCH --array=1-1648%5
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mafft_set200_%A_%a.out
#SBATCH --error=log/mafft_set200_%A_%a.err
#SBATCH --partition=scavenger

export PATH=$HOME/mafft-7.475-with-extensions/bin/:$PATH

module load NCBI-BLAST/2.7.1

busco=$(cat scripts/busco_ids_set200 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mafft-homologs.rb -l -d analyses/phylogenetics/set200/nostocales_busco_blastdb \
 -w -o '--dash --globalpair --maxiterate 100 --thread 1 --originalseqonly' \
 analyses/phylogenetics/set200/seqs/${busco}.faa > \
 analyses/phylogenetics/set200/alignments/single/${busco}_aln.faa
 