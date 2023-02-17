#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/get_trnl_set103_%A_%a.err
#SBATCH --output=log/get_trnl_set103_%A_%a.out
#SBATCH --partition=scavenger

module load NCBI-BLAST/2.7.1
module load R/4.1.1-rhel8
export PATH=/hpc/group/bio1/carlos/apps/:$PATH

bin=$(cat scripts/genome_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)
fna_file=analyses/cyano_genomes_annotation/set103/${bin}/${bin}.fna
out_dir=analyses/cyano_genomes_annotation/set103/${bin}


# Make blast database with genome contigs
makeblastdb -in ${fna_file} -input_type fasta -dbtype nucl -parse_seqids \
-out  ${out_dir}/${fna_file##*/} 
# Search trnl intron
blastn -db ${out_dir}/${fna_file##*/} -query scripts/trnl.fas \
-outfmt '6 sseqid sstart send' -out ${out_dir}/trnl_tab -max_target_seqs 1
# If there is a match...
if [[ -s ${out_dir}/trnl_tab ]] ; then
 #
 scripts/seq_range_from_blast.R ${out_dir}/trnl_tab \
 ${out_dir}/trnl_range \
 ${out_dir}/trnl_header
 # Define header and range commands for seqkit
 header=$(cat ${out_dir}/trnl_header)
 range=$(cat ${out_dir}/trnl_range)
 # Extract trnl sequence
 seqkit grep -p ${header} ${fna_file} --immediate-output |
 seqkit subseq -r ${range} > ${out_dir}/trnl.fas
 # Replace fasta header by bin id
 seq_id=$(grep ">" ${out_dir}/trnl.fas | sed 's/>//')
 sed "s/${seq_id}/${bin}/" -i ${out_dir}/trnl.fas
else
 echo ${bin}
fi
