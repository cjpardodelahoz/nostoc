
# Computational workflow for *Nostoc* project

<!-- badges: start -->
<!-- badges: end -->

This is the documentation for the data and scripts used in my paper on *Nostoc* phylogenetics, evolution and species boundaries. The manuscript is currently in prep.

You can access the tutorial for placing unknown *Nostoc* sequences within our phylogenomic framework [here](this is being built...).

## Table of Contents

## 1. How to use this guide

## 2. Computational and software setup

The analyses for this project were run on three different machines:

- [**Duke Computer Cluster**](https://oit-rc.pages.oit.duke.edu/rcsupportdocs/dcc/) (~85% of tasks): Runs on Red Hat Enterprise Linux 8 and uses the SLURM scheduler to distribute tasks.
- **iMac Pro** (~5% of tasks, mostly DiscoVista): Runs macOS V XX, Radeon Pro XX processor, 120 GB RAM.
- **MacBook Pro** (~10% of tasks, mostly R): Runs macOS Ventura 13.6.3. Apple M1, 8 Gb RAM.

I am working on containerizing all the software I used. For now, you can install things individually following the instructions in the next sections.

### 2.1 Conda environments (including R)

I exported all conda environments as YAML so you can use the exact same version and dependencies I ran. You can install each environment individually:

- [BUSCO](https://busco.ezlab.org/)
    ```sh
    # Create and install BUSCO environment
    conda env create -f conda_env_yamls/busco.yml
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
    ```
- [Quast](https://github.com/ablab/quast)
    ```sh
    conda env create -f conda_env_yamls/quast.yml
    ```
- [Graphbin2](https://github.com/metagentools/GraphBin2)
    ```sh
    # Create conda environment with dependencies
    conda env create -f conda_env_yamls/graphbin2.yml
    # Download Graphbin2 release source
    wget https://github.com/metagentools/GraphBin2/releases/download/v1.1/GraphBin2-v1.1.zip
    unzip GraphBin2-v1.1.zip
    # Add to path. This line will be in the scripts that run GraphBin2
    export PATH=$(pwd)/GraphBin2-v1.1:${PATH}
    ```
- [Prokka](https://github.com/tseemann/prokka)
    ```sh
    conda env create -f conda_env_yamls/prokka.yml
    ```
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
    ```sh
    conda env create -f conda_env_yamls/mmseqs.yml
    ```
- [Newick utilities](https://github.com/tjunier/newick_utils)
    ```sh
    conda env create -f conda_env_yamls/newick.yml
    ```
- [PlasX]()
    ```sh
    # Create conda environment with dependencies
    conda env create -f conda_env_yamls/plasx.yml
    # Close PlasX repository-There are no releases, I cloned it at commit 5c9acfd
    git clone https://github.com/michaelkyu/PlasX.git
    # Install using pip
    conda activate plasX
    cd PlasX
    pip install .
    cd ../
    # Setup PlasX
    export PATH=$(pwd)/Plasx:${PATH}
    plasx setup \
        --de-novo-families 'https://zenodo.org/record/5819401/files/PlasX_mmseqs_profiles.tar.gz' \
        --coefficients 'https://zenodo.org/record/5819401/files/PlasX_coefficients_and_gene_enrichments.txt.gz'
    ```
- [PopCOGenT](https://github.com/philarevalo/PopCOGenT)
    ```sh
    # Create conda environment with (most) dependencies
    conda env create -f conda_env_yamls/popcopgent.yml
    # Clone PopCOGenT repo–There are no releases, I cloned it at commit 7296af9
    git clone https://github.com/philarevalo/PopCOGenT.git
    cd PopCOGenT/
    conda activate PopCOGenT
    # Install required version of Infomap, distributed with popcogent
    cd Infomap/
    make
    # Download mugsy
    cd ../
    wget https://sourceforge.net/projects/mugsy/files/mugsy_x86-64-v1r2.3.tgz/download
    mv download mugsy_x86-64-v1r2.3.tgz
    tar xvzf mugsy_x86-64-v1r2.3.tgz
    cd ../ 
    ```
- [Anvio](https://anvio.org/install/)
    ```sh
    # Create conda environemtn with dependencies
    conda env create -f conda_env_yamls/anvio.yml
    # Install anvio using pip
    conda activate anvio-7.1
    curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz \
        --output ./anvio-7.1.tar.gz
    pip install ./anvio-7.1.tar.gz
    # After this is done, you will likely get an error that looks like [issue #1830](https://github.com/merenlab/anvio/issues/1839). The code below should fix it
    cp $(conda info --base)/pkgs/python-3.6.13-h12debd9_1/lib/python3.6/_sysconfigdata_x86_64_conda_linux_gnu.py \
     $(conda info --base)/envs/anvio-7.1/lib/python3.6/_sysconfigdata_x86_64_conda_linux_gnu.py
    ```
    I always ran anvio from our cluster, but I set up a port to get the interactive sessions on Chrome in my local computer. [Check out this post that explains how to do it](https://merenlab.org/2018/03/07/working-with-remote-interative/).

### 2.2 Software with pre-compiled binaries

These software have precompiled binaries. I used all of these programs on the cluster, but most of them have binaries available for macOS as well.

- [SRA toolkit](https://hpc.nih.gov/apps/sratoolkit.html)
    ```sh
    # Download toolkit
    wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    # Configure
    tar -zxf sratoolkit.3.0.0-centos_linux64.tar.gz
    cd sratoolkit.3.0.0-centos_linux64
    bin/vdb-config --interactive # Follow instructions in https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration
    ```
- [SPAdes-contig-to-graph](https://github.com/rrwick/SPAdes-Contig-Graph.git)
    ```sh
    git clone https://github.com/rrwick/SPAdes-Contig-Graph.git
    ```
- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
    ```sh
    sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
    ```
- [Seqkit](https://bioinf.shenwei.me/seqkit/)
    ```sh
    wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
    ```
- [SPAdes 3.14.1](https://github.com/ablab/spades)
    ```sh
    curl https://github.com/ablab/spades/releases/download/v3.14.1/SPAdes-3.14.1-Linux.tar.gz
    tar -zxf SPAdes-3.14.1-Linux.tar.gz
    ```
- [Bowtie2 2.3.5.1](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    ```sh
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
    ```

### 2.3 Software to build from source

These programs need to be compiled from source.

- [wASTRAL-hybrid](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral-hybrid.md)
    ```sh
    wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
    unzip Linux.zip
    cd ASTER-Linux/
    make & cd ../
    ```
- [PAML 4.9 (with MCMCTree)](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf)
    ```sh
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
    ```
- [MetaBAT 2](https://bitbucket.org/berkeleylab/metabat/src/master/)
    ```sh
    # Stable release version
    wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
    tar xzvf master.tar.gz
    cd berkeleylab-metabat-*
    #run the installation script
    mkdir build && cd build && cmake .. [ -DCMAKE_INSTALL_PREFIX=/path/to/install ] && make && make test && make install
    ```

### 2.4 Docker containers

The only program that I used as a container was [DiscoVista](https://github.com/esayyari/DiscoVista). I was not able to get this running on the cluster, so I installed it on the MacBook Pro and the iMac Pro. It worked well on both. You will need to [install Docker](https://docs.docker.com/engine/install/) for this to run. Then do:

```sh
docker pull esayyari/discovista
```

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

## 3. From raw sequence data to curated *Nostoc* genomes



#### Taxa sets

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



## 4. Phylogenetic analyses (Figure 1 and related supplements)

## 5. Evaluation of genomic species boundaries (Figure 2 and related supplements)

## 6. Analyses of *Nostoc rbcLX* sequences (Figure 3, Table 1 and related supplements)

