#!/bin/bash

app=/hpc/group/bio1/carlos/apps/PopCOGenT/src/PopCOGenT/
configfile=scripts/popcogent_config_set12c.sh
source ${configfile}
#source activate PopCOGenT
source ${mugsy_env}

if [ "${slurm_str}" = "" ]
	then
		python ${app}/get_alignment_and_length_bias.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}
		python ${app}/cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ${final_output_dir} --infomap_path ${infomap_path} ${single_cell}
fi