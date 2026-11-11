
mf_following_msas=$mf_following_msas

# variables to adapt: used for alignment
run_massivefold=$run_massivefold
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
data_dir=$data_dir
logs_dir=$logs_dir
max_template_date=$max_template_date

json_path=${input_dir}/alphafold3_json_requests/$${sequence_name}.json


date
echo "$${run_massivefold}
  --norun_inference
  --json_path $${json_path}
  --max_template_date $${max_template_date}
  --db_dir $${data_dir}
  --model_dir $${data_dir}
  --output_dir $${output}/$${sequence_name}
  --force_output_dir
"

time ${engine_prefix} $${run_massivefold} \
  --norun_inference \
  --json_path $${json_path} \
  --max_template_date $${max_template_date} \
  --db_dir $${data_dir} \
  --model_dir $${data_dir} \
  --output_dir $${output}/$${sequence_name} \
  --force_output_dir \
  || { echo "Alignment has failed. Exiting."; exit 1; }

date


if $${mf_following_msas}; then
  echo "Inference will follow after this MSAs step. Preparing batches..."
  mkdir $${output}/$${sequence_name}/$${run_name}
  json_params=$json_params
  json_msas=$${output}/$${sequence_name}/msas_alphafold3/msas_alphafold3_data.json
  batches_file=$${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
  echo "unifier
    --conversion input_inference
    --batches_file $$batches_file
    --json_params $$json_params
    --to_convert $$json_msas
    --tool AlphaFold3"

  ${mf_python} -m massivefold.parallelization.unifier \
    --conversion input_inference \
    --batches_file $$batches_file \
    --json_params $$json_params \
    --to_convert $$json_msas \
    --tool AlphaFold3 \
    || { echo "Input preparation for inference has failed. Exiting."; exit 1; }
else
  echo "Only MSAS are computed, no following inference."
fi
