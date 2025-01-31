sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
mf_plots_path=$run_massivefold_plots
scripts_dir=$scripts_dir
pkl_format=$pkl_format

# rename and move files
echo "time $${scripts_dir}/organize_outputs.py --batches_path $$output"
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

if [[ $$pkl_format == "light" ]]; then
  echo "Reducing pkl size"
  echo "$${scripts_dir}/lighten_output.py $${output}"
  time $${scripts_dir}/lighten_output.py $${output}
  echo "rm $${output}/*.pkl"
  rm $${output}/*.pkl
elif [[ $$pkl_format == "none" ]]; then
  echo "Not keeping pkl files"
  rm $${output}/*.pkl
elif [[ $$pkl_format == "full" ]]; then
  echo "Keeping pkl files as in AlphaFold2 output"
else
  echo "$$pkl_format is an invalid value for pkl_format, doing nothing"
fi

date
