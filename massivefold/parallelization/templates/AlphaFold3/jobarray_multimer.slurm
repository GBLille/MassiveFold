
# These environment variables are suggested by google deepmind at:
# https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#additional-flags
# work around a known XLA issue causing the compilation time to greatly increase
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
# token < 5,120 on aingle A100 (80 GB) or H100 (80 GB)
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
# token > 5120
: '
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=1.5
#export XLA_CLIENT_MEM_FRACTION=3.2
'
sequence_name=$sequence_name
run_name=$run_name
num_diffusion_samples=$num_diffusion_samples
data_dir=$data_dir
scripts_dir=$scripts_dir
logs_dir=$logs_dir
max_template_date=$max_template_date

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
  if [[ $$batch == 0 ]]; then
    echo "Preparing the input batches for inference..."
    json_params=$json_params
    batches_file=$${logs_dir}/$${sequence_name}/$${run_name}/$${sequence_name}_$${run_name}_batches.json
    echo "$${scripts_dir}/unifier.py
      --conversion input_inference
      --to_convert $${output_dir}/$${sequence_name}/msas_alphafold3/msas_alphafold3_data.json
      --json_params $$json_params
      --batches_file $$batches_file
      --tool AlphaFold3"

    $${scripts_dir}/unifier.py \
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
    while [[ ! -f $$path_to_msas  ]];
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
  echo "Batches had already been prepared."
fi

source $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-alphafold-3.0.1

date
echo "run_alphafold.py
  --norun_data_pipeline
  --max_template_date $${max_template_date}
  --json_path $${path_to_msas}
  --db_dir $${data_dir}
  --model_dir $${data_dir}
  --output_dir $${output_path}
  --num_diffusion_samples $${num_diffusion_samples}
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
  --jackhmmer_binary_path $$(which jackhmmer) \
  --hmmalign_binary_path $$(which hmmalign) \
  --hmmbuild_binary_path $$(which hmmbuild) \
  --hmmsearch_binary_path $$(which hmmsearch)
date
