
sequence_name=$sequence_name
run_name=$run_name
num_diffusion_samples=$num_diffusion_samples
data_dir=$data_dir
logs_dir=$logs_dir
max_template_date=$max_template_date
xla_python_client_preallocate=$xla_python_client_preallocate
xla_client_mem_fraction=$xla_client_mem_fraction
tf_force_unified_memory=$tf_force_unified_memory

batch=$$SLURM_ARRAY_TASK_ID
mf_before_inference=$mf_before_inference

json_batch_input=af3_batch_$${batch}.json
output_dir=$$(realpath ${output_dir})
output_path="$${output_dir}/$${sequence_name}/$${run_name}/"
echo output_path is $$output_path
mkdir -vp $${output_path}/
path_to_msas=$${output_dir}/$${sequence_name}/$${run_name}/$${json_batch_input}

if $${mf_before_inference}; then
  echo "Waiting for input batches for inference..."
  if [ "$$batch" = 0 ]; then
    echo "Preparing the input batches for inference..."
    json_params=$json_params
    batches_file=$${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
    echo "unifier
      --conversion input_inference
      --to_convert $${output_dir}/$${sequence_name}/msas_alphafold3/msas_alphafold3_data.json
      --json_params $$json_params
      --batches_file $$batches_file
      --tool AlphaFold3"

    unifier \
      --conversion input_inference \
      --to_convert $${output_dir}/$${sequence_name}/msas_alphafold3/msas_alphafold3_data.json \
      --json_params $$json_params \
      --batches_file $$batches_file \
      --tool AlphaFold3 \
      || { echo "Input preparation for inference has failed. Exiting."; exit 1; }
  else
    # wait until all af3_batch_n.json are finished or crash if it takes too long (10 minutes)
    timeout=600
    elapsed=0
    while [ ! -f "$$path_to_msas" ];
    do
      if (( elapsed >= timeout )); then
        echo "Something may have happened during input preparation in batch 0 (see jobarray_0.log)"
        echo "Exiting."
        exit 1
      fi
      echo "Waiting...($${elapsed}s)"
      sleep 30
      (( elapsed += 30 ))
    done
  fi
else
  echo "Batches have already been prepared."
fi

. $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-alphafold3-1.1.1

date
echo "run_alphafold.py
  --norun_data_pipeline
  --max_template_date $${max_template_date}
  --json_path $${path_to_msas}
  --db_dir $${data_dir}
  --model_dir $${data_dir}
  --output_dir $${output_path}
  --num_diffusion_samples $${num_diffusion_samples}
  --xla_python_client_preallocate $${xla_python_client_preallocate}
  --xla_client_mem_fraction $${xla_client_mem_fraction}
  --tf_force_unified_memory $${tf_force_unified_memory}
  --jackhmmer_binary_path $$(which jackhmmer)
  --hmmalign_binary_path $$(which hmmalign)
  --hmmbuild_binary_path $$(which hmmbuild)
  --hmmsearch_binary_path $$(which hmmsearch)"

time run_alphafold.py \
  --norun_data_pipeline \
  --max_template_date $${max_template_date} \
  --json_path $${path_to_msas} \
  --db_dir $${data_dir} \
  --model_dir $${data_dir} \
  --output_dir $${output_path} \
  --num_diffusion_samples $${num_diffusion_samples} \
  --xla_python_client_preallocate $${xla_python_client_preallocate} \
  --xla_client_mem_fraction $${xla_client_mem_fraction} \
  --tf_force_unified_memory $${tf_force_unified_memory} \
  --jackhmmer_binary_path $$(which jackhmmer) \
  --hmmalign_binary_path $$(which hmmalign) \
  --hmmbuild_binary_path $$(which hmmbuild) \
  --hmmsearch_binary_path $$(which hmmsearch)
date
