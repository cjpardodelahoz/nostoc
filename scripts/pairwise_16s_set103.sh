#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/pairwise_16s_set103.out
#SBATCH --error=log/pairwise_16s_set103.err
#SBATCH --partition=scavenger

# Load blast suite module
module load NCBI-BLAST/2.7.1

# Make blast database with 16s seqs
makeblastdb -in analyses/species_delimitation/16s/16s.fas \
 -input_type fasta -dbtype nucl -parse_seqids \
 -out analyses/species_delimitation/16s/16s_db
# Run all-by-all blast
blastn -query analyses/species_delimitation/16s/16s.fas \
 -db analyses/species_delimitation/16s/16s_db \
 -outfmt '6 qseqid sseqid pident' \
 -out analyses/species_delimitation/16s/blast_pairs_16s.txt
# Run all-by-all blast and print regular blast output
blastn -query analyses/species_delimitation/16s/16s.fas \
 -db analyses/species_delimitation/16s/16s_db \
 -out analyses/species_delimitation/16s/blast_pairs_16s_regfmt.txt
