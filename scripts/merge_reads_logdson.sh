#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_logdson_%A_%a.err
#SBATCH --output=log/merge_reads_logdson_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_logdson)

# Make sample specific directory
for sample in ${S} ; do
 mkdir analyses/reads/${sample}
done
# Sort reads into folders
cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib1_1_Asex_GCCAAT_L003_R1_001.fastq.gz > \
analyses/reads/JL31/JL31_R1_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib1_1_Asex_GCCAAT_L003_R2_001.fastq.gz > \
analyses/reads/JL31/JL31_R2_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib2_3_Sex_CTTGTA_L003_R1_001.fastq.gz > \
analyses/reads/JL23/JL23_R1_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib2_3_Sex_CTTGTA_L003_R2_001.fastq.gz > \
analyses/reads/JL23/JL23_R2_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib3_4_Sex_GCCAAT_L004_R1_001.fastq.gz > \
analyses/reads/JL33/JL33_R1_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib3_4_Sex_GCCAAT_L004_R2_001.fastq.gz > \
analyses/reads/JL33/JL33_R2_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib4_5_Asex_CTTGTA_L004_R1_001.fastq.gz > \
analyses/reads/JL34/JL34_R1_all.fastq.gz

cat data/reads/nic_reads/peltigera_Logdson_reads/JL3_Lib4_5_Asex_CTTGTA_L004_R2_001.fastq.gz > \
analyses/reads/JL34/JL34_R2_all.fastq.gz

