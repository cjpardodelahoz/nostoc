#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_NMS_%A_%a.err
#SBATCH --output=log/merge_reads_NMS_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_NMS)

# Make sample specific directory
for sample in ${S} ; do
 mkdir analyses/reads/${sample}
done
# Sort reads into folders
cat data/reads/nic_reads/Reads_NMS/NM-S1_S32_L008_R1_001.fastq.gz > \
analyses/reads/NMS1/NMS1_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S1_S32_L008_R2_001.fastq.gz > \
analyses/reads/NMS1/NMS1_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S2_S33_L008_R1_001.fastq.gz > \
analyses/reads/NMS2/NMS2_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S2_S33_L008_R2_001.fastq.gz > \
analyses/reads/NMS2/NMS2_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S3_S34_L008_R1_001.fastq.gz > \
analyses/reads/NMS3/NMS3_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S3_S34_L008_R2_001.fastq.gz > \
analyses/reads/NMS3/NMS3_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S4_S35_L008_R1_001.fastq.gz > \
analyses/reads/NMS4/NMS4_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S4_S35_L008_R2_001.fastq.gz > \
analyses/reads/NMS4/NMS4_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S5_S36_L008_R1_001.fastq.gz > \
analyses/reads/NMS5/NMS5_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S5_S36_L008_R2_001.fastq.gz > \
analyses/reads/NMS5/NMS5_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S6_S37_L008_R1_001.fastq.gz > \
analyses/reads/NMS6/NMS6_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S6_S37_L008_R2_001.fastq.gz > \
analyses/reads/NMS6/NMS6_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S7_S38_L008_R1_001.fastq.gz > \
analyses/reads/NMS7/NMS7_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S7_S38_L008_R2_001.fastq.gz > \
analyses/reads/NMS7/NMS7_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S8_S39_L008_R1_001.fastq.gz > \
analyses/reads/NMS8/NMS8_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S8_S39_L008_R2_001.fastq.gz > \
analyses/reads/NMS8/NMS8_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S9_S40_L008_R1_001.fastq.gz > \
analyses/reads/NMS9/NMS9_R1_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/NM-S9_S40_L008_R2_001.fastq.gz > \
analyses/reads/NMS9/NMS9_R2_all.fastq.gz

cat data/reads/nic_reads/Reads_NMS/Nostoc-R1.fastq > \
analyses/reads/NOS/NOS_R1_all.fastq

gzip analyses/reads/NOS/NOS_R1_all.fastq

cat data/reads/nic_reads/Reads_NMS/Nostoc-R2.fastq > \
analyses/reads/NOS/NOS_R2_all.fastq

gzip analyses/reads/NOS/NOS_R2_all.fastq