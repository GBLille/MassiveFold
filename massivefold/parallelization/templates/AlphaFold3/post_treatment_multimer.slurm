sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
mf_plots_path=$run_massivefold_plots
scripts_dir=$scripts_dir
batches_file=${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
pkl_format=$pkl_format

date
# adapt colabfold output for massivefold
echo "$${scripts_dir}/unifier.py
  --batches_file $${batches_file}
  --conversion output
  --to_convert $$output
  --tool AlphaFold3"

time $${scripts_dir}/unifier.py \
  --batches_file $${batches_file} \
  --conversion output \
  --to_convert $$output \
  --tool AlphaFold3

# rename and move files
echo "$${scripts_dir}/organize_outputs.py --batches_path $$output"
time $${scripts_dir}/organize_outputs.py --batches_path $$output

# plots
echo "Generating plots"

echo "python3 $$(which $${mf_plots_path})
  --input_path=$${output} 
  --top_n_predictions=$MF_plots_top_n_predictions
  --chosen_plots=$MF_plots_chosen_plots"

time python3 $$(which $${mf_plots_path}) \
  --input_path=$${output} \
  --top_n_predictions=$MF_plots_top_n_predictions \
  --chosen_plots=$MF_plots_chosen_plots \

date
