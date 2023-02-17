#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 32 # number of threads per process
#SBATCH --output=log/plasx_detect_set103_%A_%a.out
#SBATCH --error=log/plasx_detect_set103_%A_%a.err
#SBATCH --partition=scavenger

# Plasx conda env and path
source $(conda info --base)/etc/profile.d/conda.sh
conda activate plasx
export PATH=/hpc/group/bio1/carlos/apps/Plasx:${PATH}

# Genome ID variable
S=$(cat scripts/genome_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Output directory
PREFIX="analyses/plasmid_detection/set103/${S}"
# Path to genomes
GENOME_PATH="analyses/cyano_genomes/set103"
# Threads
THREADS=32
# ~1 hr if THREADS=4. ~5 min if THREADS=128.
# - For faster processing, set THREADS to the number of CPU cores available.
# - This command requires a high amount of RAM (at least ~60Gb). If your machine has low RAM, then you need to set the `--splits` flag to a high number.
#   This will split the de novo families into chunks, reducing the max RAM usage. E.g. if you have only ~8Gb, we recommend setting `--splits` to 32 or higher.
plasx search_de_novo_families \
    -g ${PREFIX}/gene_calls.txt \
    -o ${PREFIX}/de_novo_families.txt \
    --threads ${THREADS} \
    --splits 0 \
    --overwrite
# Predict plasmids
plasx predict \
    -a ${PREFIX}/cogs_and_pfams.txt ${PREFIX}/de_novo_families.txt \
    -g ${PREFIX}/gene_calls.txt \
    -o ${PREFIX}/scores.txt \
    --overwrite