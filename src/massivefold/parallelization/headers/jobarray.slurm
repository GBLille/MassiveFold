#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --array=0-$substitute_batch_number
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log

