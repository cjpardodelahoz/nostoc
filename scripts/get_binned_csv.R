#!/bin/bash

# Display help menu
if [ ${1} == "-h"] ; then
 echo "This is a script to generate comma-delimited file of bin-labeled contigs to use as input for GraphBin (bins from MetaBat2)"
 echo "Usage: `basename $0` path/to/bin/dir sample_id out_dir"
 echo "example: `basename $0` analyses/bins/metabat/8274 8274 analyses/bins/"
 exit 0
fi

bin_dir=${1}
sample_id=${2}

# Define bin number index varriable
n=1
#
while [ TEST -e ${bin_dir}/${sample_id}_bin.${n}.fa ] ; do
 grep -e ">" | sed 's/>//' > out_n
 sed -i 's/$/,'"${n}"'/' out_n
 ++n
done

cat out_* final 
 
  