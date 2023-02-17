#!/bin/bash

#SBATCH --array=1
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/get_trnl_setx_%A_%a.err
#SBATCH --output=log/get_trnl_setx_%A_%a.out
#SBATCH --partition=scavenger

module load NCBI-BLAST/2.7.1
export PATH=/hpc/group/bio1/carlos/apps/:$PATH

#bin=$(cat scripts/genome_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)
bin=3_bin.5.fa
fna_file=analyses/cyano_genomes/set8/${bin}_annotation/${bin}.fna
out_dir=analyses/cyano_genomes/set8/${bin}_annotation


# Make blast database with genome contigs
makeblastdb -in ${fna_file} -input_type fasta -dbtype nucl -parse_seqids \
-out  ${out_dir}/${fna_file##*/} 
# Search trnl intron
blastn -db ${out_dir}/${fna_file##*/} -query scripts/trnl.fas \
-outfmt '6 sseqid sseq' -out ${out_dir}/trnl_tab -max_target_seqs 1
# Convert blast output to fasta
seqkit tab2fx ${out_dir}/trnl_tab > ${out_dir}/trnl.fas
# Replace fasta header by bin id
seq_id=$(grep ">" ${out_dir}/trnl.fas | sed 's/>//')
sed "s/${seq_id}/${bin}/" -i ${out_dir}/trnl.fas
