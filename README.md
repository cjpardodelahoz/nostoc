
# *Nostoc* project data and code

<!-- badges: start -->
<!-- badges: end -->

This is the documentation for the data and scripts used in the paper by
Pardo-De la Hoz et al .

The file *command.md* all the sequence of scripts of execution for all analyses
performed in the project

## Required software

You will need to install the following ..

- R
The cluster is running on R 4.1.1-rhel8
The Local machines ran on R 4.0.5 until I started the species delimitation, then I had to upgrade to R 4.2.2 to be able to install treedataverse

```
# For R 4.2.2
# Increase memory limit for RStudio
touch .Renviron
printf "R_MAX_VSIZE = 120GB" > .Renviron
#
install.packages("tidyverse")
install.packages("rlang")
install.packages("data.table")
install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install("YuLab-SMU/treedataverse")
install.packages("labdsv")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
install.packages("ggimage")
install.packages("pegas")
install.packages("MSCquartets")
```

-conda environments

```
# BUSCO
mamba create -n busco -c bioconda -c conda-forge busco
# Download BUSCO lineage datasets for cyanobacateria and nostocales
mkdir busco_downloads
mkdir busco_downloads/lineages
wget https://busco-data.ezlab.org/v5/data/lineages/cyanobacteria_odb10.2021-02-23.tar.gz 
wget https://busco-data.ezlab.org/v5/data/lineages/nostocales_odb10.2020-03-06.tar.gz
tar -xzf nostocales_odb10.2020-03-06.tar.gz
tar -xzf cyanobacteria_odb10.2021-02-23.tar.gz
rm nostocales_odb10.2020-03-06.tar.gz
rm cyanobacteria_odb10.2021-02-23.tar.gz
mv nostocales_odb10 busco_downloads/lineages/
mv cyanobacteria_odb10 busco_downloads/lineages/
mv busco_downloads scripts/

# Quast
mamba create -n quast -c bioconda -c conda-forge quast=5.0.2

# GraphBin2
git clone https://github.com/Vini2/GraphBin2.git
cd GraphBin2/
mamba env create -f environment.yml
export PATH=/hpc/group/bio1/carlos/apps/GraphBin2:${PATH}

# Prokka
mamba create -n prokka -c conda-forge -c bioconda -c defaults prokka=1.14.6

# MMseq2s
mamba create -n mmseqs2 -c conda-forge -c bioconda mmseqs2=13.45111

# Anvi'o 7.1
# see scripts/install_anvio.sh
# I had to add the missing file _sysconfigdata_x86_64_conda_linux_gnu.py
cp /hpc/group/bio1/carlos/miniconda3_pkgs/pkgs/python-3.6.13-h12debd9_1/lib/python3.6/_sysconfigdata_x86_64_conda_linux_gnu.py python3.6/_sysconfigdata_x86_64_conda_linux_gnu.py
# Set up anvio interactive to launch from cluster and run on local browser
chmod +x scripts/dcc_interactive_setup.sh
scripts/dcc_interactive_setup.sh
# login typing dcc_interactive alias
# when logged in, set port variable
export ANVIO_PORT=8090

# Newick utilities 1.6
mamba create -n newick -c bioconda newick_utils

# PlasX (instructions from https://github.com/michaelkyu/plasx)
mamba create --name plasx -y -c anaconda -c conda-forge -c bioconda --override-channels --strict-channel-priority  numpy pandas scipy scikit-learn numba python-blosc mmseqs2=10.6d92c
git clone https://github.com/michaelkyu/PlasX.git
cd PlasX
pip install .
export PATH=/hpc/group/bio1/carlos/apps/Plasx:${PATH}
plasx setup \
    --de-novo-families 'https://zenodo.org/record/5819401/files/PlasX_mmseqs_profiles.tar.gz' \
    --coefficients 'https://zenodo.org/record/5819401/files/PlasX_coefficients_and_gene_enrichments.txt.gz'

# PopCOGenT
git clone https://github.com/philarevalo/PopCOGenT.git
cd PopCOGenT/
conda config --set restore_free_channel true
mamba env create -f PopCOGenT.yml
conda activate PopCOGenT
# Install required version of Infomap, distributed with popcogent
cd Infomap/
make
# Install mugsy
cd ../../../
wget https://sourceforge.net/projects/mugsy/files/mugsy_x86-64-v1r2.3.tgz/download
mv download mugsy_x86-64-v1r2.3.tgz
tar xvzf mugsy_x86-64-v1r2.3.tgz 
cd mugsy_x86-64-v1r2.2/
# in mugsyenv.sh modify MUGSY_INSTALL

# Gubbins
mamba create -n gubbins -c bioconda -c r -c defaults -c conda-forge gubbins
```

-To build from source:

SPAdes v. 3.14.1 https://github.com/ablab/spades/releases/tag/v3.14.1
Bowtie2 v. 2.3.5.1  https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/
MetaBAT v. 2 2020-10-09 https://bitbucket.org/berkeleylab/metabat/src/master/

-Download binaries

```
# SRA toolkit (in /hpc/group/bio1/carlos/apps)
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -zxf sratoolkit.3.0.0-centos_linux64.tar.gz
cd sratoolkit.3.0.0-centos_linux64
bin/vdb-config --interactive # Follow instructions in https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration
export PATH=/hpc/group/bio1/carlos/apps/sratoolkit.3.0.0-centos_linux64/bin:$PATH
# Repository directory is /hpc/group/bio1/carlos/apps/sratoolkit.3.0.0-centos_linux64/repo

# Script to convert contigs to assembly graph from spades-Requires Python 2.7.x
git clone https://github.com/rrwick/SPAdes-Contig-Graph.git
export PATH=/hpc/group/bio1/carlos/apps/SPAdes-Contig-Graph:$PATH

# Seqkit
wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
export PATH=/hpc/group/bio1/carlos/apps/:$PATH

# ASTRAL-III
wget https://github.com/smirarab/ASTRAL/raw/MP/Astral.5.15.5.zip
unzip Astral.5.15.5.zip
rm Astral.5.15.5.zip 
module load Java/1.8.0_60
export PATH=/hpc/group/bio1/carlos/apps/Astral:$PATH

# wASTRAL-hybrid (part of ASTER)
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
cd ASTER-Linux/
make
export PATH=/hpc/group/bio1/carlos/apps/ASTER-Linux/bin:$PATH

# PAML 4.9 (with MCMCTree)
wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
tar xf paml4.9j.tgz
rm paml4.9j.tgz
cd paml4.9j/
rm bin/*.exe
cd src
make -f Makefile
ls -lF
rm *.o
mv baseml basemlg codeml pamp evolver yn00 chi2 mcmctree ../bin
export PATH=/hpc/group/bio1/carlos/apps/paml4.9j/bin:$PATH

# Quartet scores 1.0
wget https://github.com/lutteropp/QuartetScores/releases/download/v1.0.0/QuartetScores-Linux.tar.gz
tar -zxf QuartetScores-Linux.tar.gz 
rm QuartetScores-Linux.tar.gz
export PATH=/hpc/group/bio1/carlos/apps/:$PATH
```

-Docker containers
```
# DiscoVista
# Install Docker
# Pull Discovista container
docker pull esayyari/discovista
# Download and install go
https://go.dev/doc/install
# Create a go module (which will create the go.mod file)
cd ~/gophy
go mod init gophy
# Clone the github repo
git clone https://github.com/FePhyFoFum/gophy.git
# Install the missing compiler library
brew install nlopt
# Install bp and add it to path
go build gophy/bp/bp.go
```

## Taxa sets

-set0: published *Nostoc* II genomes that were used in the nostocales paper plus the bins Nic generated from his metagenomes.

-set1: set0 plus cyanobacterial bins from my metagenome assemblies (6240_6256_7535)

-set2: all cyno bins from nic's metagenomes

-set3: published *Nostoc* II genomes that were used in the nostocales paper after removing the low quality ones based on the first report. This set also includes tw additional outgroup taxa based on the relationships in the nostocales paper: Anabaena_cylindrica_PCC_7122 and Aphanizomenon_flos_aquae_NIES_81.fa

-set4: all cyano bins from my metagenomes

-set5: *Nostoc* II bins from my metagenomes after removing low quality ones based on report 1

-set6: *Nostoc* II bins from nic after removing low quality and non-publishable ones 

-set7: *Nostoc* II bins from all newly generated metagenomes in the study after removing low quality and non-publishable ones. 

-set8: *Nostoc* II bins from all newly generated metagenomes after refinement with graphbin2

-set73: *Nostoc* II bins from sets 7 and 3 for preliminary phylogenetic analyses

-set9: *Nostoc* II bins from all newly generated metagenomesexported from anvio after visual curation.

-set10: *Nostoc* II bins from all newly generated metagenomes curated and with the extra contigs from the assembly graph.

-set103: set10 plus set3. For phylogenetic analyses.

-set103a: set10 plus set3. After sorting plasmids and chromosomes

-set103c: set10 plus set3. Putative chromosomes only

-set103p: set10 plus set3. Putative plasmids only

-set200: *Nostoc* genomes from subset 0 in Pardo-De la Hoz et al, but also including JL33. Used for dating.


