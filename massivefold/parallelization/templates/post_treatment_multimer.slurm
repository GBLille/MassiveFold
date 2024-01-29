sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
mf_plots_path=$run_massivefold_plots

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
