#!/bin/bash

#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks per node
#SBATCH --cpus-per-task=10    # Number of OpenMP threads per task
#SBATCH --hint=nomultithread # Disable hyperthreading
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log
#SBATCH --time=01:59:00      # Expected runtime HH:MM:SS (max 100h)
#SBATCH --account=$jeanzay_project@cpu #

