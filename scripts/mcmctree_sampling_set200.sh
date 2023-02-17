#!/bin/bash

#SBATCH --array=0-2
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mcmctree_sampling_set200_%A_%a.out
#SBATCH --error=log/mcmctree_sampling_set200_%A_%a.err
#SBATCH --partition=scavenger

# PAML on path
export PATH=/hpc/group/bio1/carlos/apps/paml4.9j/bin:$PATH
# Define array with chain ids
chains=('c1' 'c2' 'c3')
chain=${chains[${SLURM_ARRAY_TASK_ID}]}
# divtime directory path
divtime="analyses/phylogenetics/set200/divtime"
# npart variable
npart="1"
# part directory path
dir="${divtime}/${npart}_part"
# Alignment file name
aln="concat_ng_${npart}_part.phy"
# Tree file name
tree="set200.treefile"
# Prepare control file for mcmc sampling with data
sed -i "s|seqfile\s=\s|seqfile = ..\/..\/data\/${aln}|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|treefile\s=\s|treefile = ..\/..\/data\/${tree}|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|mcmcfile\s=\s|mcmcfile = mcmc_${chain}.txt|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|outfile\s=\s|outfile = out_${chain}.txt|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|\sndata\s=\s| ndata = ${npart}|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|seqtype\s=\s2|seqtype = 2|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|usedata\s=\s2|usedata = 2|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|burnin\s=\s50000|burnin = 200000|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|sampfreq\s=\s100|sampfreq = 100|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
sed -i "s|nsample\s=\s500000|nsample = 2000000|" \
 ${dir}/mcmc/${chain}/mcmctree.ctl
# Prepare control file for mcmc sampling from prior
sed -i "s|seqfile\s=\s|seqfile = ..\/../\/data\/${aln}|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|treefile\s=\s|treefile = ..\/../\/data/${tree}|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|mcmcfile\s=\s|mcmcfile = mcmc_pr_${chain}.txt|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|outfile\s=\s|outfile = out_pr_${chain}.txt|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|\sndata\s=\s| ndata = ${npart}|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|seqtype\s=\s2|seqtype = 2|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|usedata\s=\s2|usedata = 0|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|burnin\s=\s50000|burnin = 200000|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|sampfreq\s=\s100|sampfreq = 100|" \
 ${dir}/prior/${chain}/mcmctree.ctl
sed -i "s|nsample\s=\s500000|nsample = 2000000|" \
 ${dir}/prior/${chain}/mcmctree.ctl
# Copy the file with hessian and gradient to the mcmc and prior dirs
cp ${dir}/gH/rst2 ${dir}/mcmc/${chain}/in.BV
cp ${dir}/gH/rst2 ${dir}/prior/${chain}/in.BV
# Go to mcmc directory and run mcmc sampling
cd ${dir}/mcmc/${chain}
mcmctree mcmctree.ctl
# Go to prior directory and run prior sampling
cd ../../prior/${chain}
mcmctree mcmctree.ctl
