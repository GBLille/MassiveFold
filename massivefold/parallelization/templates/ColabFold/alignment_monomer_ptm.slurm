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
  echo "ValueError: --pair_strategy\n'$$pair_strategy'"
  exit 1
fi

echo "Pairing strategy used: $$pairing_strategy"

msas_output="$${output}/$${sequence_name}/msas/"

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
