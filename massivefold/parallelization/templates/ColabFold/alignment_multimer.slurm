# variables to adapt: used for alignment
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir
scripts_dir=$scripts_dir

pair_strategy=$pair_strategy
if [[ $$pair_strategy == "greedy" ]]; then
  pairing_strategy=0
elif [[ $$pair_strategy == "complete" ]]; then
  pairing_strategy=1
else
  echo "ValueError: pair_strategy should be 'greedy' or 'complete'"
  exit 1
fi

echo "./$${scripts_dir}/unifier.py
  --conversion input
  --to_convert $$fafile"

./$${scripts_dir}/unifier.py \
  --conversion input \
  --to_convert $$fafile

fafile=${input_dir}/converted_for_colabfold/$${sequence_name}.fasta

msas_output="$${output}/$${sequence_name}/msas_colabfold/"

date

echo "colabfold_search
  $${fafile}
  $${data_dir}
  $${msas_output}
  --pairing_strategy $${pairing_strategy}"

time colabfold_search \
  $${fafile} \
  $${data_dir} \
  $${msas_output} \
  --pairing_strategy $${pairing_strategy}

date
