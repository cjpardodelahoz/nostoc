#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/rbclx_abmi_epa_placement.out
#SBATCH --error=log/rbclx_abmi_epa_placement.err
#SBATCH --partition=scavenger

# Load RAxML module
module load RAxML/8.2.12
# Path to wd. Replace with the path where you have the nostoc project dir
# This is because RAxML doesn't like relative paths
wd="/hpc/group/bio1/carlos/nostoc"
# Run placement
raxmlHPC-PTHREADS-SSE3 -f v \
 -s analyses/species_delimitation/rbclx/clade_assignment/alignments/rbclx_set103_abmi_public_aln.phy \
 -t analyses/species_delimitation/rbclx/clade_assignment/trees/placement/backbone.tree \
 -w ${wd}/analyses/species_delimitation/rbclx/clade_assignment/trees/placement \
 -m GTRCAT -n epa_result -T 4
