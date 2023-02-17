#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/plasx_anvi_setup.out
#SBATCH --error=log/plasx_anvi_setup.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

THREADS=8

# Download COGs (~2 min on fast network) and Pfam (~20 min)
# - The flags `--cog-version COG14` and `--pfam-version 32.0` directs anvio to download the 2014 version of the COG database and v32.0 of Pfam, which are used by PlasX. Without these flags, anvio will by default download the latest versions of COG (v. 2020) and Pfam (v35.0), which PlasX does not use.
# - The flags `--cog-data-dir COG2014` and `--pfam-data-dir Pfam_v32` directs anvio to store the COG and Pfam databases in new subfolders named `COG_2014` and `Pfam_v32`, respectively.
anvi-setup-ncbi-cogs --cog-version COG14 --cog-data-dir scripts/COG_2014 -T ${THREADS}
anvi-setup-pfams --pfam-version 32.0 --pfam-data-dir scripts/Pfam_v32