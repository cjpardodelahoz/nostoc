#!/bin/bash

# This is a script to create the clade_def file that includes all bipartitions in a
# reference tree. When the clade_def file is generated with this function, the bipartitions will
# be in the same order as R assigns node labels when reading the tree, which allows connecting the 
# bipartitions with the node labels. Bipartition 1 will have all the taxa in the tree. 
# Bipartition b will correspond to internal node n+b in R, where n is the number of taxa in the tree.

# Run in shell as:
# ./get_clade_def.sh path/to/input/tree out/dir outgroupfile

# Get the bipartitions from the tree using the R script
chmod +x scripts/get_biparts.R 
scripts/get_biparts.R ${1} ${2}/biparts.temp ${3}
# Reformat the bipartitions for discovista
cat ${2}/biparts.temp | \
 awk '{gsub("\\,", "\"\"\+\"\""); print}' | \
 sed 's/^/"/' > ${2}/biparts
# Generate the clade-def file for disco vista from the biparts file
touch ${2}/clade_def
printf "Clade Name\tClade Definition\tSection Letter\tComponents\tShow\tComments\n" > \
	${2}/clade_def
clade=1
for bipart in $(< ${2}/biparts) ; do
	printf "${clade}\t${bipart}\"\"\"\tNone\t\t1\t\n" >> \
	${2}/clade_def
	let "clade++"
done
# Relabel the first bipartition as the "All" row that contains the labels of all taxa
sed -i'.bak' 's/^1\t/All\t/' ${2}/clade_def
rm ${2}/clade_def.bak
