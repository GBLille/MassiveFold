#!/bin/bash

# slurm parameters to launch on cpu

# activate AF/MF conda env
# e.g: module load massivefold/1.0.0

module purge
module load massivefold/1.0.0

sequence_name=$sequence_name
run_name=$run_name
output=../output_array/$${sequence_name}/$${run_name}/
mf_plots_path= #/gpfswork/rech/fvp/commun/plots/MF_plots.py


# rename and move files
echo "time ./organize_outputs.py --batches_path $$output"
time ./organize_outputs.py --batches_path $$output

# plots
echo "Generating plots"

echo "python3 $${mf_plots_path}
        --input_path=$${output} 
        --top_n_predictions=$MF_plots_top_n_predictions
        --chosen_plots=$MF_plots_chosen_plots"

time python3 $${mf_plots_path} \
        --input_path=$${output} \
        --top_n_predictions=$MF_plots_top_n_predictions \
        --chosen_plots=$MF_plots_chosen_plots \

date
