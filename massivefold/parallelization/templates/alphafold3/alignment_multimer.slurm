
# variables to adapt: used for alignment
run_massivefold=$run_massivefold
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir
scripts_dir=$scripts_dir

json_params=./alphafold3_params.json
model_dir=../../data/af3weights 

echo "./$${scripts_dir}/unifier.py
  --conversion input
  --json_params $$json_params
  --to_convert $$fafile
  --tool alphafold3"

./$${scripts_dir}/unifier.py \
  --conversion input \
  --json_params $$json_params \
  --to_convert $$fafile \
  --tool alphafold3

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
