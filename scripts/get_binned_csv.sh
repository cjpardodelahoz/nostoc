#!/bin/bash

# Display help menu
if [ ${1} == "-h"] ; then
 echo "This is a script to generate comma-delimited file of bin-labeled contigs to use as input for GraphBin (bins from MetaBat2)"
 echo "Usage: `basename $0` path/to/bin/dir sample_id out_dir"
 echo "example: `basename $0` analyses/bins/metabat/8274 8274 analyses/bins/graphbin/binned_csv"
 exit 0
fi

# Define variables from arguments
bin_dir=${1}
sample_id=${2}
out_dir=${3}

# Make output directoey
mkdir -p ${outdir}

# Define bin number index varriable
n=1
#
while [ TEST -e ${bin_dir}/${sample_id}_bin.${n}.fa ] ; do
 grep -e ">" | sed 's/>//' > ${out_dir}/${sample_id}_bin.${n}
 sed -i 's/$/,'"${n}"'/' ${out_dir}/${sample_id}_bin.${n}
 ++n
done
#
cat ${out_dir}/${sample_id}_bin.${n} ${out_dir}/${sample_id}_contig_labels.csv 
 
  