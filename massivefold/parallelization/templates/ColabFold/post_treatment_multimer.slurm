sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
mf_plots_path=$run_massivefold_plots
scripts_dir=$scripts_dir
batches_file=${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
pkl_format=$pkl_format

# adapt colabfold output for massivefold
echo "$${scripts_dir}/unifier.py
  --batches_file $${batches_file}
  --conversion output
  --to_convert $$output
  --tool ColabFold
  "

time $${scripts_dir}/unifier.py \
  --batches_file $${batches_file} \
  --conversion output \
  --to_convert $$output \
  --tool ColabFold

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

if [[ $$pkl_format == "light" ]]; then
  pkl_mode="light"
  echo "Reducing pkl size"
elif [[ $$pkl_format == "none" ]]; then
  echo "Not keeping pkl files"
  pkl_mode="delete"
elif [[ $$pkl_format == "full" ]]; then
  echo "Keeping pkl files as in AlphaFold output"
  pkl_mode="full"
else
  echo "$$pkl_format is an invalid value for pkl_format, doing nothing"
fi

echo "$${scripts_dir}/lighten_output.py $${output} --pickle_size $${pkl_mode}"
time $${scripts_dir}/lighten_output.py $${output} --pickle_size $${pkl_mode}

date

mkdir $${output}/log
output_logs_path=$$(realpath $${output}/log)
cp -r ${logs_dir}/$${sequence_name}/$${run_name}/* $${output}/log/
cd $${output}/log 
zip ../logs.zip *
rm $$output_logs_path/*
mv ../logs.zip .
