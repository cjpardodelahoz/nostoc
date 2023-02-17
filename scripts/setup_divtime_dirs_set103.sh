#!/bin/bash

part_schemes="1_part 4_part"
mcmctree_template_path="analyses/phylogenetics/set103/mcmctree_templates"
divtime_path="analyses/phylogenetics/set103/divtime"
tree="analyses/phylogenetics/set103/trees/astral/wastral_divtime.phy"
for scheme in ${part_schemes} ; do
 mkdir -p ${divtime_path}/${scheme}/data
 mkdir -p ${divtime_path}/${scheme}/gH
 mkdir -p ${divtime_path}/${scheme}/mcmc/c{1..3}
 mkdir -p ${divtime_path}/${scheme}/prior/c{1..3}
 cp analyses/phylogenetics/set103/alignments/concat/concat_ng_${scheme}.phy \
  ${divtime_path}/${scheme}/data
 cp ${tree} ${divtime_path}/${scheme}/data
done
for scheme in ${part_schemes} ; do
# cp ${mcmctree_template_path}/lg.dat ${divtime_path}/${scheme}/gH
 cp ${mcmctree_template_path}/mcmctree-outBV.ctl ${divtime_path}/${scheme}/gH
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/mcmc/c1
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/mcmc/c2
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/mcmc/c3
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/prior/c1
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/prior/c2
 cp ${mcmctree_template_path}/mcmctree.ctl ${divtime_path}/${scheme}/prior/c3
done