#!/bin/bash

#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks per node
#SBATCH --cpus-per-task=10    # Number of OpenMP threads per task
#SBATCH --hint=nomultithread # Disable hyperthreading
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/post_treatment_%j.log
#SBATCH --time=00:05:00      # Expected runtime HH:MM:SS (max 100h)
#SBATCH --account=$jeanzay_project@cpu #

module purge
module load massivefold/1.0.0

sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
run_massivefold_plots=$run_massivefold_plots

# rename and move files
echo "time ./organize_outputs.py --batches_path $$output"
time ./organize_outputs.py --batches_path $$output

# plots
echo "Generating plots"

echo "python3 $$run_massivefold_plots   
        --input_path=$${output} 
        --top_n_predictions=$MF_plots_top_n_predictions
        --chosen_plots=$MF_plots_chosen_plots"

time python3 $$run_massivefold_plots \
        --input_path=$${output} \
        --top_n_predictions=$MF_plots_top_n_predictions \
        --chosen_plots=$MF_plots_chosen_plots \
