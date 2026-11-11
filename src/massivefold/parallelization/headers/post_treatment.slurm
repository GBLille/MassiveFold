#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log

