sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
mf_plots_path=$run_massivefold_plots
scripts_dir=$scripts_dir
keep_pkl=$keep_pkl
batches_file=${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json

# adapt colabfold output for massivefold
echo "$${scripts_dir}/colabfold_adapter.py
  --batches_file $${batches_file}
  --conversion output
  --to_convert $$output"

time $${scripts_dir}/colabfold_adapter.py \
  --batches_file $${batches_file} \
  --conversion output \
  --to_convert $$output

# rename and move files
echo "$${scripts_dir}/organize_outputs.py --batches_path $$output"
time $${scripts_dir}/organize_outputs.py --batches_path $$output

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

if ! $$keep_pkl; then
  rm $${output}/*pkl
fi

date
