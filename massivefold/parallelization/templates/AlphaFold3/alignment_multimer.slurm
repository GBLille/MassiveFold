
mf_following_msas="$${MF_FOLLOWING_MSAS}"

# variables to adapt: used for alignment
run_massivefold=$run_massivefold
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
data_dir=$data_dir
logs_dir=$logs_dir
scripts_dir=$scripts_dir

source $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-alphafold3
json_path=${input_dir}/alphafold3_json_requests/$${sequence_name}.json

date
echo "$${run_massivefold}
  --norun_inference
  --json_path $${json_path}
  --db_dir $${data_dir}
  --model_dir $${data_dir}
  --output_dir $${output}/$${sequence_name}
  --jackhmmer_binary_path $$(which jackhmmer)
  --hmmalign_binary_path $$(which hmmalign)
  --hmmbuild_binary_path $$(which hmmbuild)
  --hmmsearch_binary_path $$(which hmmsearch)
"

time $${run_massivefold} \
  --norun_inference \
  --json_path $${json_path} \
  --db_dir $${data_dir} \
  --model_dir $${data_dir} \
  --output_dir $${output}/$${sequence_name} \
  --jackhmmer_binary_path $$(which jackhmmer) \
  --hmmalign_binary_path $$(which hmmalign) \
  --hmmbuild_binary_path $$(which hmmbuild) \
  --hmmsearch_binary_path $$(which hmmsearch)
date

conda deactivate

if [[ -n $${mf_following_msas}  ]] && [[ $${mf_following_msas} == "true"  ]]; then
  echo "Inference will folow after this MSAS step. Preparing batches..."
  mkdir $${output}/$${sequence_name}/$${run_name}
  json_params=$json_params
  json_msas=$${output}/$${sequence_name}/msas_alphafold3/msas_alphafold3_data.json
  batches_file=$${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
  echo "./$${scripts_dir}/unifier.py
    --conversion input_inference
    --batches_file $$batches_file
    --json_params $$json_params
    --to_convert $$json_msas
    --tool AlphaFold3"

  ./$${scripts_dir}/unifier.py \
    --conversion input_inference \
    --batches_file $$batches_file \
    --json_params $$json_params \
    --to_convert $$json_msas \
    --tool AlphaFold3
else
  echo "Only MSAS are computed, no following inference."
fi
