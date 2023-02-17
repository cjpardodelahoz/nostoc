#!/bin/bash

#SBATCH --array=1-121
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/copy_curated_to_set9_%A_%a.out
#SBATCH --error=log/copy_curated_to_set9_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

# Variable with Sample id
S=$(cat scripts/sample_ids_anvio_edited | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variable with bin ids
bins=$(ls analyses/bins/anvio/${S}/summary_out/bin_by_bin)

# Copy bins to the set9 directory
# I set this up so that the names of the bin files are built from the sample IDS
# and using only the bin number from the anvio export. This is to minimize potential
# mistakes I may have made when manualy entering the bin names in anvi-refine.
for bin in ${bins} ; do
 cp analyses/bins/anvio/${S}/summary_out/bin_by_bin/${bin}/${bin}-contigs.fa \
 analyses/cyano_genomes/set9/${S}_bin_${bin##*_}.fa
done
