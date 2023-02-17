#!/bin/bash

#SBATCH --array=1,4
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/hessian_set103_%A_%a.out
#SBATCH --error=log/hessian_set103_%A_%a.err
#SBATCH --partition=scavenger

# PAML on path
export PATH=/hpc/group/bio1/carlos/apps/paml4.9j/bin:$PATH
# divtime directory path
divtime="analyses/phylogenetics/set103/divtime"
# npart variable
npart="${SLURM_ARRAY_TASK_ID}"
# part directory path
dir="${divtime}/${npart}_part"
# Alignment file name
aln="concat_ng_${npart}_part.phy"
# Tree file name
tree="wastral_divtime.phy"
# Prepare control file for Hessian
sed -i "s|seqfile\s=\s|seqfile = ..\/data\/${aln}|" \
 ${dir}/gH/mcmctree-outBV.ctl
sed -i "s|treefile\s=\s|treefile = ..\/data\/${tree}|" \
 ${dir}/gH/mcmctree-outBV.ctl
sed -i "s|\sndata\s=\s| ndata = ${npart}|" \
 ${dir}/gH/mcmctree-outBV.ctl
sed -i "s|seqtype\s=\s|seqtype = 0|" \
 ${dir}/gH/mcmctree-outBV.ctl
# Get control file for gradient and Hessian estimate
cd ${dir}/gH # Has to run from that directory to write the outputs there
mcmctree mcmctree-outBV.ctl
