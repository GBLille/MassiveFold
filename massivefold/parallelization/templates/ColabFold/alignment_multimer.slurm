# variables to adapt: used for alignment
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir
scripts_dir=$scripts_dir

pair_strategy=$pair_strategy
templates=$templates
custom_template_path=$custom_template_path

echo "./$${scripts_dir}/unifier.py
  --conversion input
  --to_convert $$fafile"

./$${scripts_dir}/unifier.py \
  --conversion input \
  --to_convert $$fafile

fafile=${input_dir}/converted_for_colabfold/$${sequence_name}.fasta

msas_output="$${output}/$${sequence_name}/msas/"

date

echo "colabfold_search
  $${fafile}
  $${data_dir}
  $${msas_output}"

time colabfold_search \
  $${fafile} \
  $${data_dir} \
  $${msas_output}

date
