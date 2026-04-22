
 pdb_db_templates=pdb100_230517
# variables to adapt: used for alignment
use_templates=$templates
sequence_name=$sequence_name
run_name=$run_name
output=$output_dir
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir

pair_strategy=$pair_strategy
if [ "$$pair_strategy" = "greedy" ]; then
  pairing_strategy=0
elif [ "$$pair_strategy" = "complete" ]; then
  pairing_strategy=1
else
  echo "ValueError: --pair_strategy\n'$$pair_strategy'"
  exit 1
fi

echo "Pairing strategy used: $$pairing_strategy"

msas_output="$${output}/$${sequence_name}/msas_colabfold/"

date

. $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-colabfold-1.6.1

echo "colabfold_search
  $${fafile}
  $${data_dir}
  $${msas_output}
  --pairing_strategy $${pairing_strategy}
  --use-templates $${use_template}
  --db2 $${pdb_db_templates}"


time colabfold_search \
  $${fafile} \
  $${data_dir} \
  $${msas_output} \
  --pairing_strategy $${pairing_strategy} \
  --use-templates $${use_template} \
  --db2 $${pdb_db_templates}

date
