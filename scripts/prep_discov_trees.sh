#!/bin/bash

# This is a script to prepare trees in the input format for Discovista. Run in shell as:
# prep_discov_trees.sh path/to/trees extension path/to/discov_in type

for tree in $(cd ${1} && ls *${2}) ; do
  mkdir -p ${3}/${tree}-${4}
  cp ${1}/${tree} ${3}/${tree}-${4}/estimated_species_tree.tree
done