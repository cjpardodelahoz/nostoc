#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/pairwise_all_rbclx.out
#SBATCH --error=log/pairwise_all_rbclx.err
#SBATCH --partition=scavenger

# Load blast suite module
module load NCBI-BLAST/2.7.1

# MAke directory to store blast results
mkdir -p analyses/species_delimitation/rbclx/clade_assignment/blast
# Make blast database with 16s seqs
makeblastdb -in analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_set103_abmi_public.fna \
 -input_type fasta -dbtype nucl -parse_seqids \
 -out analyses/species_delimitation/rbclx/clade_assignment/blast/rbclx_db
# Run all-by-all blast
blastn -query analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_set103_abmi_public.fna \
 -db analyses/species_delimitation/rbclx/clade_assignment/blast/rbclx_db \
 -outfmt '6 qseqid sseqid pident length gaps' \
 -out analyses/species_delimitation/rbclx/clade_assignment/blast/blast_pairs_rbclx.txt
# Run all-by-all blast and print regular blast output
blastn -query analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_set103_abmi_public.fna \
 -db analyses/species_delimitation/rbclx/clade_assignment/blast/rbclx_db \
 -out analyses/species_delimitation/rbclx/clade_assignment/blast/blast_pairs_rbclx_regfmt.txt
