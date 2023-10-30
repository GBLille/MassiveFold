#!/bin/bash

# slurm parameters to launch on cpu

# activate AF/MF conda env
# e.g: module load massivefold/1.0.0

mf_plots_path=
sequence_name=$jobname
output=../output_array

# rename and move files
echo "time ./organize_outputs.py --batches_path $$output/$$sequence_name"
time ./organize_outputs.py --batches_path $$output/$$sequence_name

# delete not needed files
rm -r $$output/$$sequence_name/batch_*

# plots
echo "Generating plots"
echo "python3 $${mf_plots_path}
  --input_path $${output}/$${sequence_name}
  --plot_type for_each
  --top_n_predictions 5
  --chosen_plots coverage,DM_plddt_PAE"


time python3 $${mf_plots_path} \
  --input_path $${output}/$${sequence_name} \
  --plot_type for_each \
  --top_n_predictions 5  \
  --chosen_plots coverage,DM_plddt_PAE

date
