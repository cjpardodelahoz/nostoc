#!/bin/bash

################################################################################
# ASSEMBLE AND BIN METAGENOMES
################################################################################

# Assemble Novaseq order 7535

# make copy of reads
cp data/reads/Carlos_7535_220221A6/*.gz analyses/reads/
# Generate file with sample ids
files=$(cd analyses/reads && ls *.gz)
for file in ${files} ; do
 echo ${file%_S*} >> scripts/files
done
cat scripts/files | uniq > scripts/sample_ids
rm scripts/files
# Organize reads
gunzip analyses/reads/*gz # make into a slurm script
sbatch scripts/merge_reads.sh
# Get QC reports with fastqc
sbatch scripts/fastqc_array.sh
# Trim reads with Trimmomatic
sbatch scripts/trimmomatic_array.sh
# Assemble metagenomes
sbatch scripts/spades_array.sh 

# Assemble Novaseq orders 6240 and 6526

# Generate files with sample ids
s6240=$(cd data/reads/6240_1_2/ && ls *.gz)
s6526=$(cd data/reads/6526 && ls *.gz)
for file in ${s6240} ; do
 echo ${file%_S*} >> scripts/s6240
done
cat scripts/s6240 | uniq > scripts/sample_ids_6240
for file in ${s6526} ; do
 echo ${file%_S*} >> scripts/s6526
done
cat scripts/s6526 | uniq > scripts/sample_ids_6526
cat scripts/s6240 scripts/s6526 | sort | uniq > scripts/sample_ids_6240_6526
rm scripts/s6240
rm scripts/s6526
# Copy and merge and sort the 6240 reads
sbatch scripts/merge_reads_6240.sh 
# Copy and sort 6526 reads
sbatch scripts/merge_reads_6526.sh 
# Get QC reports with fastqc
sbatch scripts/fastqc_6240_6526.sh
# Trim reads with Trimmomatic
sbatch scripts/trimmomatic_6240_6526.sh
# Get file with all sample ids
cat scripts/sample_ids_6240_6526 scripts/sample_ids > scripts/sample_ids_6240_6526_7535
# Get post-trim QC report of all reads
sbatch scripts/fastqc_trim_6240_6526_7535.sh
# Remove raw reads from analyses directory
sbatch scripts/delete_untrimmed_reads.sh
# Assemble metagenomes
sbatch scripts/spades_6240_6526.sh 

# Assemble nic reads [assemblies in work since 9/6/2022]

# Generate files with sample ids (note that the sample ids of resequencing 2021
# are a subset of the ids of genomes 2020)
genomes2020=$(cd data/reads/nic_reads/genomes2020 && ls *.gz)  # genomes 2020 and resequening 2021
for file in ${genomes2020} ; do
 echo ${file%_S*} | sed 's/6126-//' >> scripts/genomes2020_resequencing2021
done
cat scripts/genomes2020_resequencing2021 | sort| uniq > \
scripts/sample_ids_genomes2020_resequencing2021
rm scripts/genomes2020_resequencing2021
intermediaries=$(cd data/reads/nic_reads/reads_genomes_intermediaires/ && ls *.gz) # genomes from intermediaries
for file in ${intermediaries} ; do
 echo ${file%_S*} >> scripts/intermediaries
done
cat scripts/intermediaries | sort | uniq > scripts/sample_ids_intermediares
# Copy and merge all nic reads
sbatch scripts/merge_reads_genomes2020_resequencing2021.sh
sbatch scripts/merge_reads_intermediaries.sh
sbatch scripts/merge_reads_logdson.sh
sbatch scripts/merge_reads_NMS.sh
# Get file with sample ids for all of nic samples
cat scripts/sample_ids_genomes2020_resequencing2021 \
scripts/sample_ids_intermediares \
scripts/sample_ids_logdson \
scripts/sample_ids_NMS > \
scripts/sample_ids_nic
# Get QC reports with fastqc
sbatch scripts/fastqc_nic.sh
# Trim reads with Trimmomatic
sbatch scripts/trimmomatic_nic.sh
# Get post-trim QC report of all reads
sbatch scripts/fastqc_trim_nic.sh
# Remove raw reads from analyses directory
sbatch scripts/delete_untrimmed_reads.sh
# Assemble metagenomes
sbatch scripts/spades_nic.sh 

# Get and assemble reads from Vecherskii et al. 2022

# Set path for SRA toolkits
export PATH=/hpc/group/bio1/carlos/apps/sratoolkit.3.0.0-centos_linux64/bin:$PATH
# Download reads from SRA
prefetch -O ./ $(< scripts/vecherskii_accesions) 
fasterq-dump -O ./ $(< scripts/vecherskii_accesions)#
# Move reads to data directory
mkdir data/reads/vecherskii
mv SR* data/reads/vecherskii
# Copy reads to analyses directory
sbatch scripts/copy_vecherskii.sh
# Get QC reports with fastqc
sbatch scripts/fastqc_vecherskii.sh
# Trim reads with Trimmomatic
sbatch scripts/trimmomatic_vecherskii.sh
# Get post-trim QC report of all reads
sbatch scripts/fastqc_trim_vecherskii.sh
# Remove raw reads from analyses directory
sbatch scripts/delete_untrimmed_reads.sh
# Assemble metagenomes
sbatch scripts/spades_vecherskii.sh 

# Bin Novaseq assemblies 6240, 6526 and 7535

# Tar and gzip copy of reads
sbatch scripts/targz_6240_6526_7535.sh 
# Build assembly database
sbatch scripts/bt2_db_6240_6526_7535.sh 
# Map assembly to corrected reads
sbatch scripts/bt2_map_6240_6526_7535.sh
# Convert SAM file to sorted BAM
sbatch scripts/samtosortedbam_6240_6526_7535.sh
# Remove SAM file to save storage space
sbatch scripts/remove_sam_6240_6526_7535.sh 
# Summarize BAM depth
sbatch scripts/sumdepth_6240_6526_7535.sh
# Bin contigs with MetaBat
mkdir analyses/bins
mkdir analyses/bins/metabat
for S in $(cat scripts/sample_ids_6240_6526_7535) ; do
 mkdir analyses/bins/metabat/${S}
done
sbatch scripts/metabat_6240_6526_7535.sh

# Bin assemblies from nic

# Build assembly database
sbatch scripts/bt2_db_nic.sh 
# Map assembly to corrected reads
sbatch scripts/bt2_map_nic.sh
# Convert SAM file to sorted BAM
sbatch scripts/samtosortedbam_nic.sh
# Remove SAM file to save storage space
sbatch scripts/remove_sam_nic.sh 
# Summarize BAM depth
sbatch scripts/sumdepth_nic.sh
# Bin contigs with MetaBat
sbatch scripts/metabat_nic.sh

# Bin assemblies from vecherskii

# Build assembly database
sbatch scripts/bt2_db_vecherskii.sh 
# Map assembly to corrected reads
sbatch scripts/bt2_map_vecherskii.sh
# Convert SAM file to sorted BAM
sbatch scripts/samtosortedbam_vecherskii.sh
# Remove SAM file to save storage space
sbatch scripts/remove_sam_vecherskii.sh 
# Summarize BAM depth
sbatch scripts/sumdepth_vecherskii.sh
# Bin contigs with MetaBat
for S in $(cat scripts/sample_ids_vecherskii) ; do
 mkdir analyses/bins/metabat/${S}
done
sbatch scripts/metabat_vecherskii.sh

################################################################################
# TAXON SAMPLING
################################################################################

# Get bin taxonomy with CheckM-6240, 6526 and 7535
sbatch scripts/checkm_6240_6526_7535.sh
# Get paths to cyanobacteria bins
chmod 777 scripts/get_cyano_bin_paths_6240_6256_7535.R
scripts/get_cyano_bin_paths_6240_6256_7535.R 
# Make copy of cyano bins from 6240_6256_7535 into set1 directory
metabat_bin_paths=$(cat scripts/bin_paths_6240_6256_7535)
for bin in ${metabat_bin_paths} ; do
 cp ${bin} analyses/cyano_genomes/set1/${bin##*/}
done
# Copy genomes from set0 into set1 directory and add .fa extension
set0=$(ls analyses/cyano_genomes/set0/ | sed 's/.fna//')
for bin in ${set0} ; do
 cp analyses/cyano_genomes/set0/${bin}* analyses/cyano_genomes/set1/${bin}.fa
done
# Run GUNC on set1 genomes
sbatch scripts/gunc_set1.sh
# Run checkm on set1 genomes
sbatch scripts/checkm_set1.sh
# Run BUSCO with cyano_odb_10 on set1 genomes
sbatch scripts/busco_set1.sh
# Sort BUSCO seqs
scripts/sort_busco_seqs_set1.R
# Align amino acids for set1 single genes
sbatch scripts/mafft_set1_aa.sh 
# Align nucleotides for set1 using aa alignments as guide
sbatch scripts/pal2nal_set1.sh
# Trim gaps from alignments
sbatch scripts/trimal_ng_set1.sh
# remove trimal crap from headers
sbatch scripts/fix_headers_set1.sh
# Concatenate sequences
sbatch scripts/concatenate_set1_na.sh
# Run concatenated trees
sbatch scripts/concat_tree_set1_na.sh
sbatch scripts/concat_tree_set1_aa.sh
# Get report number 1
scripts/report_1.R

# Let's bring in nic's reassembled genomes

# Get bin taxonomy with CheckM
sbatch scripts/checkm_nic.sh
# Get paths to cyanobacteria bins
chmod +x scripts/get_cyano_bin_paths_nic.R
scripts/get_cyano_bin_paths_nic.R 
# Make copy of cyano bins from nic into set2 directory
mkdir -p analyses/cyano_genomes/set2
metabat_bin_paths=$(cat scripts/bin_paths_nic)
for bin in ${metabat_bin_paths} ; do
 cp ${bin} analyses/cyano_genomes/set2/${bin##*/}
done

# Sort bins into new taxa sets

# Make taxon set3
cp analyses/cyano_genomes/set1/* analyses/cyano_genomes/set3/
rm analyses/cyano_genomes/set3/P*
rm analyses/cyano_genomes/set3/S*
rm analyses/cyano_genomes/set3/NMS*
rm analyses/cyano_genomes/set3/JL*
rm analyses/cyano_genomes/set3/Nostoc_sp_UIC_10630.fa
rm analyses/cyano_genomes/set3/X2*
rm analyses/cyano_genomes/set3/3_*
rm analyses/cyano_genomes/set3/8274*
rm analyses/cyano_genomes/set3/8277*
# Make taxon set 4
mkdir analyses/cyano_genomes/set4
metabat_bin_paths=$(cat scripts/bin_paths_6240_6256_7535)
for bin in ${metabat_bin_paths} ; do
 cp ${bin} analyses/cyano_genomes/set4/${bin##*/}
done
# Make taxon set5
mkdir analyses/cyano_genomes/set5
metabat_bin_paths=$(cat scripts/kept_bin_paths_6240_6256_7535)
for bin in ${metabat_bin_paths} ; do
 cp ${bin} analyses/cyano_genomes/set5/${bin##*/}
done
# Make taxon set 6
mkdir analyses/cyano_genomes/set6
metabat_bin_paths=$(cat scripts/kept_bin_paths_nic)
for bin in ${metabat_bin_paths} ; do
 cp ${bin} analyses/cyano_genomes/set6/${bin##*/}
done
# Make taxon set 7
mkdir analyses/cyano_genomes/set7
cp analyses/cyano_genomes/set5/* analyses/cyano_genomes/set7
cp analyses/cyano_genomes/set6/* analyses/cyano_genomes/set7

# Curation of newly-generated cyano genomes (The Making of set9 and set10)

# genome id's for set7 (129 with multi strain bins)
ls analyses/cyano_genomes/set7 > scripts/genome_ids_set7
# Get sample (metagenome) ids for bins in set 7 NO DUPLICATES HERE (127)
cat scripts/genome_ids_set7 | sed 's/_bin.*//' | sort | uniq > scripts/sample_ids_set7
# Run quast and busco for set 7
sbatch scripts/quast_set7.sh
sbatch scripts/busco_set7.sh
# Get contigs assembly graph from spades output
sbatch scripts/contigs_to_graph_set7.sh
# Label contigs for graphbin
sbatch scripts/label_contigs_set7.sh
# Refine bins with grahbin
sbatch scripts/graphbin_set7.sh
sbatch scripts/graphbin_to_collection_set7.sh
# Bin contigs from graphbin csv output # module load R/4.1.1-rhel8
chmod +x scripts/bin_from_csv_set7.R
sbatch scripts/bin_from_csv_set7.sh
# Make taxon set 8
mkdir analyses/cyano_genomes/set8
cp scripts/genome_ids_set7 scripts/genome_ids_set8
genome_ids_set7=$(cat scripts/genome_ids_set7)
for bin in ${genome_ids_set7} ; do
 S=$(echo ${bin} | sed 's/_.*//')
 cp analyses/bins/graphbin/${S}/${bin} analyses/cyano_genomes/set8/${bin}
done
# Run quast gunc busco and CheckM for set 8
sbatch scripts/quast_set8.sh
sbatch scripts/checkm_set8.sh
sbatch scripts/gunc_set8.sh
sbatch scripts/busco_set8.sh
# Get qc table for set 8
chmod +x scripts/report_2.R
sbatch scripts/report_2.sh
chmod +x scripts/compare_busco_set7_set8.R
scripts/compare_busco_set7_set8.R
# Annotate set 8 genomes with Prokka
sbatch scripts/prokka_set8.sh
# Get 16s from annotated genomes set8 FIX THIS TO ACCOUNT FOR MULTIPLE AND PARTIAL HITS
sbatch scripts/get_16s_from_ffn_set8.sh
# MAKE PYTHON TOOL TO DETECT tRNAs in PROKKA OUTPUT
# trnL gene is an un-annotated copy of the tRNA-Leu(taa)!
# Get trnl from annotated assemblies for set8
chmod +x scripts/seq_range_from_blast.R
sbatch scripts/get_trnl_set8.sh
# Get csv files of removed and added contigs after graphbin-in bandage csv format
sbatch scripts/get_delta_contigs_set8.sh
# Get anvio collection with graphbin delta
sbatch scripts/graphbin_to_collection_set8.sh
# set up anvio contigs database for qc
sbatch scripts/anvio_contigs_db_set8.sh
sbatch scripts/index_bam_set8.sh
sbatch scripts/anvio_profile_set8.sh
sbatch scripts/anvio_import_graphbin_collection_set8.sh
# Annotate contig taxonomy with mmseqs2 (UniRef90 database)
sbatch scripts/download_uniref90.sh
sbatch scripts/mmseqs_query_dbs_set8.sh
sbatch scripts/mmseqs_taxonomy_set8.sh
sbatch scripts/mmseqs_taxonomy_tsv_set8.sh
# Create file with taxonomy and delta-contigs info to visualize with anvio
chmod 777 scripts/get_delta_contigs_taxonomy_layer.R
sbatch scripts/get_delta_contigs_taxonomy_layer_set8.sh
# Launch anvi-refine
# I edited graphbin_detal_collection.txt and delta_contigs_taxonomy_layer.tsv 
# manually to remove dulpicated contig names that were assigned to two bins.
# when logged in, set port variable
export ANVIO_PORT=8090
conda activate anvio-7.1
S=P8570
anvi-refine --profile-db analyses/bins/anvio/${S}/profile_1/PROFILE.db \
 --contigs-db analyses/bins/anvio/${S}/contigs.db \
 --collection-name graphbin_delta \
 --bin-ids-file analyses/bins/anvio/${S}/bin_ids_graphbin_delta \
 --additional-layers analyses/bins/anvio/${S}/delta_contigs_taxonomy_layer.tsv
# Export bins edited with anvio with anvi-summarize
sbatch scripts/anvio_export_curated_cyano_bins.sh
# Copy curated bins into set9 directory
mkdir -p analyses/cyano_genomes/set9
sbatch scripts/copy_curated_to_set9.sh
# Add Nostoc bin from P9820 which was mislabeled in anvio
# cp analyses/bins/anvio/P9820/summary_out/bin_by_bin/P8920_bin_6/P8920_bin_6-contigs.fa \
# analyses/cyano_genomes/set9/P9820_bin_6.fa
# Add the right bin for P10247
cp analyses/cyano_genomes/set4/P10247_bin.16.fa analyses/cyano_genomes/set9/P10247_bin_16.fa
# Make set10 by adding the extra contigs from the assembly graph
mkdir -p analyses/cyano_genomes/set10
cp analyses/cyano_genomes/set9/* analyses/cyano_genomes/set10
sbatch scripts/add_extra_contigs_from_graph_set9.sh
# Add the bin from P8690 which comes entirely from the graph to set10
cp analyses/genome_qc/set8/delta_contigs/P8690_bin.1.fa/extra_contigs_from_graph.fasta \
 analyses/cyano_genomes/set10/P8690_bin_1.fa
# Remove extra characters from bandage from P8690_bin_1.fa
sed -i 's/+//' analyses/cyano_genomes/set10/P8690_bin_1.fa
sed -i 's/-//' analyses/cyano_genomes/set10/P8690_bin_1.fa


################################################################################
# PHYLOGENETIC INFERENCE
################################################################################

# Make taxon set 103
mkdir -p analyses/cyano_genomes/set103
cp analyses/cyano_genomes/set10/* analyses/cyano_genomes/set103
cp analyses/cyano_genomes/set3/* analyses/cyano_genomes/set103
rm analyses/cyano_genomes/set103/P8840_bin_nostoc.fa 
ls analyses/cyano_genomes/set103 > scripts/genome_ids_set103
# Run BUSCO with nostocales_odb_10 on set73 genomes
sbatch scripts/busco_set103.sh
# Sort BUSCO seqs
chmod +x scripts/sort_busco_seqs_set103.R
sbatch scripts/sort_busco_seqs_set103.sh
# Annotate set 103 genomes with Prokka
sbatch scripts/prokka_set103.sh
# Get 16s from annotated genomes set8 FIX THIS TO ACCOUNT FOR MULTIPLE AND PARTIAL HITS
sbatch scripts/get_16s_from_ffn_set103.sh
# Get trnl from annotated assemblies for set103
chmod +x scripts/seq_range_from_blast.R
sbatch scripts/get_trnl_set103.sh
# Compile 16s and trnl seqs and align manually. 16s from Nostoc_KVJ2.fa was
# removed because it wasn't nostoc
sbatch scripts/compile_trnl_16s_set103.sh
# Align amino acids for set103 single genes using mafft-dash
# busco locus 4681at1161 was excluded from downstreama analyses because it contained
# sequences that appear fragmented and no sites were left after removing gaps
mkdir -p analyses/phylogenetics/set103/alignments/single
sbatch scripts/mafft_set103_aa.sh 
# Align nucleotides for set1 using aa alignments as guide
sbatch scripts/pal2nal_set103.sh
# Trim gaps from alignments
sbatch scripts/trimal_ng_set103.sh
# Remove trimal crap from headers
sbatch scripts/fix_headers_set103.sh
# Summarize single gene alignments after trimming
mkdir -p analyses/phylogenetics/set103/alignments/sumaries
sbatch scripts/summarize_na_alignments_set103.sh
# Filter busco ids to remove loci with < 200 variable sites and < 136 (90%) taxa
# This resulted in 1517 loci being kept
chmod +x scripts/filter_busco_ids_set103.R
sbatch scripts/filter_busco_ids_set103.sh
# Get codon partitions
sbatch scripts/get_codon_partition_set103.sh
# Run gene trees
mkdir -p analyses/phylogenetics/set103/trees/single
sbatch scripts/ml_gene_trees_set103.sh
sbatch scripts/ml_16s_trnl_trees_set103.sh
# Replace empty ufboot values for 0s in the single gene trees
chmod +x scripts/replace_empty_ufboot_single_set103.R
Rscript scripts/replace_empty_ufboot_single_set103.R
# Get Majority Rule Consensus tree from the set103 gene trees
chmod +x scripts/get_mr_tree_set103.R
Rscript scripts/get_mr_tree_set103.R
# Run ASTRAL on all gene trees
sbatch scripts/astral_set103.sh
# Run ASTRAL on gene trees with branches < %10 UFBoot collapsed
sbatch scripts/astral10_set103.sh
# Run ASTRAL weighted by UFBoot support and branch lengths
sbatch scripts/wastral_set103.sh
# Concatenate sequences
sbatch scripts/concatenate_set103.sh
# Get codon partition file for the concatenated alignment
# Remember to remove the codon lines for 16s and trnl
sbatch scripts/get_codon_partition_concat_set103.sh
# Run concatenated trees
sbatch scripts/concat_pf_set103.sh
sbatch scripts/concat_tree_set103.sh

################################################################################
# PHYLOGENETIC CONFLICT ANALYSES
################################################################################

# Gene trees vs Astral

# LOCALLY ON TITUS
# Prep gene trees for discovista
scripts/prep_discov_trees.sh analyses/phylogenetics/set103/trees/single \
 noempty.treefile analyses/phylogenetics/set103/conflict/single_vs_astral/discovista_in single
# Get clade definition file for Discovista using the concat tree as reference
chmod +x scripts/get_clade_def.sh
scripts/get_clade_def.sh analyses/phylogenetics/set103/trees/astral/astral.tree \
 analyses/phylogenetics/set103/conflict/single_vs_astral \
 scripts/outgroup_set103
# Run Discovista
cp -r analyses/phylogenetics/set103/conflict/single_vs_astral ~/
cd ~
docker run -v /Users/lutzonilab/single_vs_astral10:/data esayyari/discovista discoVista.py -m 0 -k 1 \
 -c ./clade_def \
 -p ./discovista_in/ \
 -t 95 -o ./discovista_out
 
# Gene trees vs Astral10 (nodes < 10% UFBoot collapsed)

# LOCALLY ON FRANK
# Prep gene trees for discovista
scripts/prep_discov_trees.sh analyses/phylogenetics/set103/trees/single \
 noempty.treefile analyses/phylogenetics/set103/conflict/single_vs_astral10/discovista_in single
# Get clade definition file for Discovista using the concat tree as reference
chmod +x scripts/get_clade_def.sh
scripts/get_clade_def.sh analyses/phylogenetics/set103/trees/astral/astral10.tree \
 analyses/phylogenetics/set103/conflict/single_vs_astral10 \
 scripts/outgroup_set103
# Run Discovista
cp -r analyses/phylogenetics/set103/conflict/single_vs_astral10 ~/
cd ~
docker run -v /Users/lutzonilab/single_vs_astral10:/data esayyari/discovista discoVista.py  -m 0 -k 1 \
 -c ./clade_def \
 -p ./discovista_in/ \
 -t 95 -o ./discovista_out
 
# Gene trees vs wAstral-hybrid

# LOCALLY ON FRANK
# Prep gene trees for discovista
scripts/prep_discov_trees.sh analyses/phylogenetics/set103/trees/single \
 noempty.treefile analyses/phylogenetics/set103/conflict/single_vs_wastral/discovista_in single
# Get clade definition file for Discovista using the concat tree as reference
chmod +x scripts/get_clade_def.sh
scripts/get_clade_def.sh analyses/phylogenetics/set103/trees/astral/wastral.tree \
 analyses/phylogenetics/set103/conflict/single_vs_wastral \
 scripts/outgroup_set103
# Run Discovista
cp -r analyses/phylogenetics/set103/conflict/single_vs_wastral ~/
cd ~
docker run -v /Users/lutzonilab/single_vs_wastral:/data esayyari/discovista discoVista.py -m 0 -k 1 \
 -c ./clade_def \
 -p ./discovista_in/ \
 -t 95 -o ./discovista_out
 
# Gene trees vs Concat

# LOCALLY ON FRANK
# Prep gene trees for discovista
scripts/prep_discov_trees.sh analyses/phylogenetics/set103/trees/single \
 noempty.treefile analyses/phylogenetics/set103/conflict/single_vs_concat/discovista_in single
# Get clade definition file for Discovista using the concat tree as reference
chmod +x scripts/get_clade_def.sh
scripts/get_clade_def.sh analyses/phylogenetics/set103/trees/concat/concat_tree_ng_na.treefile \
 analyses/phylogenetics/set103/conflict/single_vs_concat \
 scripts/outgroup_set103
# Run Discovista
cp -r analyses/phylogenetics/set103/conflict/single_vs_concat ~/
cd ~
docker run -v /Users/lutzonilab/single_vs_concat:/data esayyari/discovista discoVista.py -m 0 -k 1 \
 -c ./clade_def \
 -p ./discovista_in/ \
 -t 95 -o ./discovista_out
 
# Quartet scores
 
# Compute quartet scores using the weighted astral tree as a reference
# and raw as well as collapsed (<95%) gene trees
sbatch scripts/quartet_scores_set103.sh

# Preliminary phylogenetic analyses of conflict

# Make taxon set 73
mkdir -p analyses/cyano_genomes/set73
cp analyses/cyano_genomes/set7/* analyses/cyano_genomes/set73
cp analyses/cyano_genomes/set3/* analyses/cyano_genomes/set73
ls analyses/cyano_genomes/set73 > scripts/genome_ids_set73
# Run BUSCO with cyano_odb_10 on set73 genomes
sbatch scripts/busco_set73.sh
# Sort BUSCO seqs
chmod +x scripts/sort_busco_seqs_set73.R
sbatch scripts/sort_busco_seqs_set73.sh
# Align amino acids for set73 single genes
mkdir -p analyses/phylogenetics/set73/alignments/single
sbatch scripts/mafft_set73_aa.sh 
# Align nucleotides for set1 using aa alignments as guide
sbatch scripts/pal2nal_set73.sh
# Trim gaps from alignments
sbatch scripts/trimal_ng_set73.sh
# Remove trimal crap from headers
sbatch scripts/fix_headers_set73.sh
# Run gene trees
mkdir -p analyses/phylogenetics/set73/trees/single/
sbatch scripts/ml_gene_trees_set73.sh
# Concatenate sequences
sbatch scripts/concatenate_set73.sh
# Run concatenated trees
sbatch scripts/concat_tree_set73.sh
# LOCALLY ON TITUS
# Prep gene trees for discovista
scripts/prep_discov_trees.sh analyses/phylogenetics/set73/trees/single \
treefile analyses/phylogenetics/set73/conflict/discovista_in single
# Get clade definition file for Discovista using the concat tree as reference
# This still needed some manual modifications in the clade_def file
# Will figure it out later
# Missing 155	"3_bin.5.fa""+""Nostoc_sp_Peltigera_membranacea_cyanobiont_210A.fa"""	None		1	
chmod +x scripts/get_clade_def.sh
scripts/get_clade_def.sh analyses/phylogenetics/set73/trees/concat/concat_ng_na.treefile \
analyses/phylogenetics/set73/conflict
# Run Discovista
cp -r analyses/phylogenetics/set73/conflict ~/
cd ~
for tree in $(ls conflict/discovista_in) ; do
 sed -i'.bak' 's/):/)0:/g' conflict/discovista_in/${tree}/estimated_species_tree.tree
 rm conflict/discovista_in/${tree}/*.bak
done
docker run -v /Users/titus/conflict:/data esayyari/discovista discoVista.py -m 0 -k 1 \
 -c ./clade_def \
 -p ./discovista_in/ \
 -t 95 -o ./discovista_out
# MOVE conflict directory back
# Get conflict pies
scripts/conflict_pies_set73.R

################################################################################
# DIVERGENCE TIME ESTIMATION
################################################################################

# Time estimates of the early divergences of Nostoc in te context of nostocales

# Directory for set200 aa sequences (1648 busco loci from Pardo-De la Hoz et al)
mkdir -p analyses/phylogenetics/set200/seqs
# Extract set200 seqs 
sbatch scripts/extract_seqs_set200.sh
# Align amino acid sequences
mkdir -p analyses/phylogenetics/set200/alignments/single
sbatch scripts/mafft_set200_aa.sh 
# Trim gaps from alignments
sbatch scripts/trimal_ng_set200.sh
# Remove trimal crap from headers
sbatch scripts/fix_headers_set200.sh
# Concatenate sequences and generate 3 partition schemes:
# Single partition with all genes (concat_ng_1_part.phy) 
# 1648 partitions with each gene as a separate one (concat_ng_1648_part.phy)
# 4 partitions each of 412 randome genes (concat_ng_4_part.phy)
sbatch scripts/concatenate_set200.sh
# Directories for divergence time estimation
chmod +x scripts/setup_divtime_dirs.sh
scripts/setup_divtime_dirs.sh
# Get gradient, Hessian and branch lenghts for mcmc sampling (for 1 partition only)
sbatch scripts/hessian_set200.sh
# Run mcmc sampling from prior and posterior
sbatch scripts/mcmctree_sampling_set200.sh

# Divergence time estimation in Nostoc

# Concatenate sequences and generate 2 partition schemes:
# Single partition with all genes (concat_ng_1_part.phy) 
# 4 partitions each of 412 randome genes (concat_ng_4_part.phy)
sbatch scripts/concatenate_divtime_set103.sh
# Prepare wASTRAL tree for divergence time estimation by removing labels and 
# branch lengths
Rscript scripts/prep_wastral_divtime_set103.R
# Manually add thecalibration points using the results from set200
# Directories for divergence time estimation
chmod +x scripts/setup_divtime_dirs_set103.sh
scripts/setup_divtime_dirs_set103.sh
# Get gradient, Hessian and branch lenghts for mcmc sampling (for 1 partition only)
sbatch scripts/hessian_set103.sh
# Run mcmc sampling from prior and posterior
sbatch scripts/mcmctree_sampling_1_part_set103.sh
sbatch scripts/mcmctree_sampling_4_part_set103.sh

# Plots for Figure 1 - time tree and conflict sumaries
Rscript scripts/fig1.R

################################################################################
# PLASMID DETECTION
################################################################################

# Dowlonad the pfam and COG databases that anvio needs for gene annotation
sbatch scripts/plasx_anvi_setup.sh
# Annotate COGs and Pfams and export to text files using anvio. This will also
# simplify the fasta headers and write the contigs to contigs.fasta
sbatch scripts/plasx_anvi_annotate_set103.sh
# Predict plasmids with plasx
sbatch scripts/plasx_detect_set103.sh
# Predict plasmids with deeplasmid (locally on FRANK)
scripts/deeplasmid_set103.sh
# Compile depth files to use for contig classification
mkdir analyses/plasmid_detection/set103/depths
for genome in $(cat misc_files/genome_ids_set103) ; do 
 cp analyses/assemblies/${genome%_bin_*}/${genome%_bin_*}_assembly_depths.txt \
  analyses/plasmid_detection/set103/depths/${genome}_depths.txt
done
# Summarize plasmid detection results and sort contigs by class for each genome
sort_plasmids_set103.R
# Sort genomes into three new sets: 
# set103a (all 103 including plasmids and chromosomes)
# set103c (chromosomes only)
# set103p (plasmids only)
mkdir -p analyses/cyano_genomes/set103a
mkdir -p analyses/cyano_genomes/set103c
mkdir -p analyses/cyano_genomes/set103p
for genome in $(cat misc_files/genome_ids_set103) ; do
 cp analyses/plasmid_detection/set103/${genome}/${genome%.fa}*.fa \
  analyses/cyano_genomes/set103a/
 cp analyses/plasmid_detection/set103/${genome}/${genome%.fa}_chromosome.fa \
  analyses/cyano_genomes/set103c/
 cp analyses/plasmid_detection/set103/${genome}/${genome%.fa}_plasmid.fa \
  analyses/cyano_genomes/set103p/
done

################################################################################
# SPECIES DELIMITATION
################################################################################

# FastANI PopCOGenT and 16S clustering

# Run all-by-all comparison of set12c with Fastani
mkdir -p analyses/species_delimitation/fastani/set12c
sbatch scripts/fastani_set12c.sh
# PopCOGenT on chromosomes to test recombination
sbatch scripts/popcogent_set12c.sh
# Run all-by-all 16S comparisons
mkdir -p analyses/species_delimitation/16s
cp analyses/phylogenetics/set103/seqs/16s.fas analyses/species_delimitation/16s/16s.fas
sbatch scripts/pairwise_16s_set103.sh
# Summarize clustering results
# This will generate the plots summarizing the distribution of pairwise ANI,
# the correlation of ANI with alignment fraction, and the ANI gap
Rscript scripts/sum_clusters.R

# Recombination with Gubbins
conda activate gubbins
generate_ska_alignment.py --reference analyses/cyano_genomes/set12c/Nostoc_sp_Peltigera_membranacea_cyanobiont_232_chromosome.fa \
 --input analyses/species_delimitation/gubbins/test/test_df1.tsv \
 --out analyses/species_delimitation/gubbins/test/out1.aln \
 --threads 4

run_gubbins.py --prefix analyses/species_delimitation/gubbins/test/out1 \
 --threads 4 \
 analyses/species_delimitation/gubbins/test/out1.aln
 
generate_ska_alignment.py --reference analyses/cyano_genomes/set12c/Nostoc_sp_Peltigera_membranacea_cyanobiont_232_chromosome.fa \
 --input analyses/species_delimitation/gubbins/test/test_df2.tsv \
 --out analyses/species_delimitation/gubbins/test/out2.aln \
 --threads 4

run_gubbins.py --prefix analyses/species_delimitation/gubbins/test/out2 \
 --threads 4 \
 analyses/species_delimitation/gubbins/test/out2.aln

# Testing differation in co-occurrence in three lineage case studies with rbcLX

# Make directory structure
mkdir -p analyses/species_delimitation/cooccurrence/seqs
mkdir -p analyses/species_delimitation/cooccurrence/alignments
mkdir -p analyses/species_delimitation/cooccurrence/trees/placement
# Copy rbcL (2400at1161) and rbcX sequences (14372at1161) from set103 taxa
# I also made a copy of the rbcLX seqs from ABMI sites under seqs/
cp analyses/phylogenetics/set103/seqs/2400at1161.fna \
 analyses/species_delimitation/cooccurrence/seqs/rbcl_set103.fna
cp analyses/phylogenetics/set103/seqs/14372at1161.fna \
 analyses/species_delimitation/cooccurrence/seqs/rbcx_set103.fna
# Concatenate rbcl and rbcx from set103
sbatch scripts/concatenate_rbclx.sh
# Add the ABMI  and global rbclx seqs to set103
# I edited the ABMI seqs to exclude the universal tags and primers and removed
# P11480, and put them all in the right orientation (rbclx_abmi_edited_all.fna)
cat analyses/species_delimitation/cooccurrence/seqs/rbclx_set103.fna \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_abmi_edited_all.fna \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_global_v.fna \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_global_xxxix.fna> \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_abmi_global.fna
# Align all the rbcLX
# I edited this alignment in Mesquite  and exluded ambiguous regions
sbatch scripts/mafft_rbclx_placement.sh
# Place the rbclx
sbatch scripts/rbclx_abmi_epa_placement.sh
# Sort rbclx seqs into focal groups
Rscript scripts/sort_rbclx_placements.R
# Align rbcLX focal groups
sbatch scripts/mafft_rbclx_focal_groups.sh
# Infer ML trees for focal groups
sbatch scripts/ml_rbclx_focal_groups.sh
# Plot trees for focal groups
Rscript scripts/plot_focal_trees.R




# 
Rscript scripts/fastbaps_cluster.R


################################################################################
# VOUCHER TABLE
################################################################################

# Get list of genomes curated and generated in this study to clean
# Remobe P8840_bin_nostoc.fa from the list
ls analyses/cyano_genomes/set10/ > misc_files/genome_ids_set10
# Run NCBI Foreign Contamination Screen pipeline to clean set103c and set103p
sbatch scripts/ncbi_fcs_set103c.sh
sbatch scripts/ncbi_fcs_set103p.sh
# Set11c with newy generated genomes, chromosomes only
mkdir analyses/cyano_genomes/set11c
for genome in $(cat misc_files/genome_ids_set10) ; do
 cp analyses/genome_qc/set103c/fcs/${genome%.fa}_chromosome.fa/cleaned_sequences/${genome%.fa}_chromosome.fa \
  analyses/cyano_genomes/set11c/${genome%.fa}_chromosome.fa
 sed -i".bak" "s|lcl\|||" analyses/cyano_genomes/set11c/${genome%.fa}_chromosome.fa
done
rm analyses/cyano_genomes/set11c/*.bak
# Set11p with newy generated genomes, plasmids only
mkdir analyses/cyano_genomes/set11p
for genome in $(cat misc_files/genome_ids_set10) ; do
 cp analyses/genome_qc/set103p/fcs/${genome%.fa}_plasmid.fa/cleaned_sequences/${genome%.fa}_plasmid.fa \
  analyses/cyano_genomes/set11p/${genome%.fa}_plasmid.fa
 sed -i".bak" "s|lcl\|||" analyses/cyano_genomes/set11p/${genome%.fa}_plasmid.fa
done
rm analyses/cyano_genomes/set11p/*.bak
# set12c set12p with set11cp and plasmid-sorted plublic genomes from set103
mkdir -p analyses/cyano_genomes/set12c
mkdir -p analyses/cyano_genomes/set12p
cp analyses/cyano_genomes/set11c/* analyses/cyano_genomes/set12c
cp analyses/cyano_genomes/set11p/* analyses/cyano_genomes/set12p
ls analyses/cyano_genomes/set3 > misc_files/genome_ids_set3
for genome in $(cat misc_files/genome_ids_set3) ; do
 cp analyses/cyano_genomes/set103c/${genome%.fa}_chromosome.fa \
  analyses/cyano_genomes/set12c
 cp analyses/cyano_genomes/set103p/${genome%.fa}_plasmid.fa \
  analyses/cyano_genomes/set12p
done
# set12c and set12p genome file list
ls analyses/cyano_genomes/set12c > misc_files/genome_ids_set12c
ls analyses/cyano_genomes/set12p > misc_files/genome_ids_set12p
# Run Busco, GUNC and Quast on set12c
sbatch scripts/busco_set12c.sh
sbatch scripts/quast_set12c.sh
sbatch scripts/gunc_set12c.sh
#
Rscript scripts/report_voucher.R

################################################################################
# GLOBAL RBCLX RECLASSIFICATION
################################################################################

# Data compilation

# Compile public rbclx metadata
Rscript scripts/compile_public_rbclx_metadata.R

# PULL
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/scripts .
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/document .
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/species_delimitation analyses
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/cyano_genomes analyses
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/genome_qc/set103 analyses/genome_qc
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/phylogenetics/set103/conflict analyses/phylogenetics/set103
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/assemblies analyses
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/cyano_genomes_annotation analyses
rsync -av cjp47@dcc-login.oit.duke.edu:/work/cjp47/nostoc/assemblies analyses



# PUSH
rsync -av scripts cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc
rsync -av packrat cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc
rsync -av .Rprofile cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/.Rprofile
rsync -av analyses/cyano_genomes cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses
rsync -av analyses/phylogenetics/set103/trees cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/phylogenetics/set103
rsync -av analyses/bins/graphbin cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc/analyses/bins
