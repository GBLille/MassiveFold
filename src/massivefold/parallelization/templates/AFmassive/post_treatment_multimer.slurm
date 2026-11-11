sequence_name=$sequence_name
run_name=$run_name
output=${output_dir}/$${sequence_name}/$${run_name}/
batches_file=${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
pkl_format=$pkl_format

# adapt colabfold output for massivefold
echo "unifier
  --batches_file $${batches_file}
  --conversion output
  --to_convert $$output
  --tool AFmassive
  "

time ${mf_python} -m massivefold.parallelization.unifier \
  --batches_file $${batches_file} \
  --conversion output \
  --to_convert $$output \
  --tool AFmassive \
  || { echo "Output standardization has failed. Exiting."; exit 1; }

# rename and move files
echo "time ${mf_python} -m massivefold.parallelization.organize_outputs --batches_path $$output --nomultithreading"
time ${mf_python} -m massivefold.parallelization.organize_outputs --batches_path $$output --nomultithreading || { echo "Batch gathering has failed. Exiting."; exit 1; }

# plots
echo "Generating plots"

echo "massivefold_plots
        --input_path=$${output} 
        --top_n_predictions=$MF_plots_top_n_predictions
        --chosen_plots=$MF_plots_chosen_plots"

time ${mf_python} -m massivefold.massivefold_plots \
        --input_path=$${output} \
        --top_n_predictions=$MF_plots_top_n_predictions \
        --chosen_plots=$MF_plots_chosen_plots \
        || echo "Plots production has failed."

if [ "$$pkl_format" = "light" ]; then
  pkl_mode="light"
  echo "Reducing pkl size"
elif [ "$$pkl_format" = "none" ]; then
  echo "Not keeping pkl files"
  pkl_mode="delete"
elif [ "$$pkl_format" = "full" ]; then
  echo "Keeping pkl files as in AlphaFold output"
  pkl_mode="full"
else
  echo "$$pkl_format is an invalid value for pkl_format, doing nothing"
fi

echo "lighten_output $${output} --pickle_size $${pkl_mode}"
time ${mf_python} -m massivefold.parallelization.lighten_output $${output} --pickle_size $${pkl_mode} || { echo "Output lightening failed. Exiting"; exit 1; }

date

mkdir $${output}/log
output_logs_path=$$(realpath $${output}/log)
cp -r ${logs_dir}/$${sequence_name}/$${run_name}/* $${output}/log/
cd $${output}/log 
${mf_python} -m zipfile -c ../logs.zip *
rm $$output_logs_path/*
mv ../logs.zip .
