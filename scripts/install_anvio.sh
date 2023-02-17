#!/bin/bash

#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/install_anvio.out
#SBATCH --error=log/install_anvio.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh

conda create -y --name anvio-7.1 python=3.6
conda activate anvio-7.1

mamba install -y -c bioconda "sqlite>=3.31.1"
mamba install -y -c bioconda prodigal
mamba install -y -c bioconda mcl
mamba install -y -c bioconda muscle=3.8.1551
mamba install -y -c bioconda hmmer
mamba install -y -c bioconda diamond
mamba install -y -c bioconda blast
mamba install -y -c bioconda megahit
mamba install -y -c bioconda spades
mamba install -y -c bioconda bowtie2 tbb=2019.8
mamba install -y -c bioconda bwa
mamba install -y -c bioconda samtools=1.9
mamba install -y -c bioconda centrifuge
mamba install -y -c bioconda trimal
mamba install -y -c bioconda iqtree
mamba install -y -c bioconda trnascan-se
mamba install -y -c bioconda r-base
mamba install -y -c bioconda r-stringi
mamba install -y -c bioconda r-tidyverse
mamba install -y -c bioconda r-magrittr
mamba install -y -c bioconda r-optparse
mamba install -y -c bioconda bioconductor-qvalue
mamba install -y -c bioconda fasttree
mamba install -y -c bioconda vmatch
mamba install -y -c bioconda fastani

curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz \
        --output ../apps/anvio-7.1.tar.gz
        
pip install ../apps/anvio-7.1.tar.gz