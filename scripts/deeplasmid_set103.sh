#!/usr/bin/bash 

# Copy plasmid detection directory to root to run deeplasmid
cp -r analyses/plasmid_detection ~/plasmid_detection
# Iterate through the genomes and run deeplasmid
for genome in $(cat scripts/genome_ids_set103 | tail -n 63) ; do
 docker run -it -v ~/plasmid_detection/set103/${genome}/contigs.fasta:/srv/jgi-ml/classifier/dl/in.fasta \
  -v ~/plasmid_detection/set103/${genome}:/srv/jgi-ml/classifier/dl/outdir \
  billandreo/deeplasmid feature_DL_plasmid_predict.sh in.fasta outdir
done

 