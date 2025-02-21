source $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-afmassive-1.1.5

date

#variables to adapt
sequence_name=$sequence_name
run_name=$run_name
fafile=${input_dir}/$${sequence_name}.fasta
run_massivefold=$run_massivefold
data_dir=$data_dir
scripts_dir=$scripts_dir
max_template_date=$max_template_date
#max_template_date=2023-07-31
uniref_database=$uniref_database

if [ -z "$$uniref_database" ]; then
  echo "No specified path for uniref30, use default ${data_dir}/uniref30/UniRef30_2021_03"
  uniref_database=$${data_dir}/uniref30/UniRef30_2021_03
fi

output_dir=$$(realpath ${output_dir})
output_path="$${output_dir}/$${sequence_name}/$${run_name}/batch_$${SLURM_ARRAY_TASK_ID}"
echo output_path is $$output_path
mkdir -vp $${output_path}/$${sequence_name}
ln -s $${output_dir}/$${sequence_name}/msas $${output_path}/$${sequence_name}

# fixed variables
use_precomputed_msas=true
alignments_only=false
use_gpu_relax=false
models_to_relax=none

# variable as MF_parallel parameters
model_preset=$model_preset
db_preset=$db_preset
dropout=$dropout
dropout_structure_module=$dropout_structure_module
dropout_rates_filename=$dropout_rates_filename
max_recycles=$max_recycles
early_stop_tolerance=$early_stop_tolerance
bfd_max_hits=$bfd_max_hits
mgnify_max_hits=$mgnify_max_hits
uniprot_max_hits=$uniprot_max_hits
uniref_max_hits=$uniref_max_hits
templates=$templates
min_score=$min_score
max_score=$max_score
stop_recycling_below=$stop_recycling_below

# parameters for each task in the job array
batches_path=${logs_dir}/${sequence_name}/${run_name}
batch_model=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element model)
batch_start=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element start)
batch_end=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element end)

echo Model used is $$batch_model, batch computed between prediction $$batch_start and $$batch_end

echo "python $$(which $${run_massivefold})
    --fasta_paths=$${fafile}
    --output_dir=$${output_path}
    --data_dir=$${data_dir}
    --db_preset=$${db_preset}
    --model_preset=$${model_preset}
    --models_to_relax=$${models_to_relax}
    --use_gpu_relax=$${use_gpu_relax}
    --max_template_date=$${max_template_date}
    --use_precomputed_msas=$${use_precomputed_msas}
    --end_prediction=$${batch_end}
    --alignments_only=$${alignments_only}
    --dropout=$${dropout}
    --dropout_structure_module=$${dropout_structure_module}
    --dropout_rates_filename=$${dropout_rates_filename}
    --max_recycles=$${max_recycles}
    --early_stop_tolerance=$${early_stop_tolerance}
    --min_score=$${min_score}
    --max_score=$${max_score}
    --stop_recycling_below=$${stop_recycling_below}
    --bfd_max_hits=$${bfd_max_hits}
    --mgnify_max_hits=$${mgnify_max_hits}
    --uniprot_max_hits=$${uniprot_max_hits}
    --uniref_max_hits=$${uniref_max_hits}
    --models_to_use=$${batch_model}
    --start_prediction=$${batch_start}
    --templates=$${templates}
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb70_database_path=$${data_dir}/pdb70/pdb70
    --uniref30_database_path=$${uniref_database}"

time python $$(which $${run_massivefold}) \
    --fasta_paths=$${fafile} \
    --output_dir=$${output_path} \
    --data_dir=$${data_dir} \
    --db_preset=$${db_preset} \
    --model_preset=$${model_preset} \
    --models_to_relax=$${models_to_relax} \
    --use_gpu_relax=$${use_gpu_relax} \
    --max_template_date=$${max_template_date} \
    --use_precomputed_msas=$${use_precomputed_msas} \
    --end_prediction=$${batch_end} \
    --alignments_only=$${alignments_only} \
    --dropout=$${dropout} \
    --dropout_structure_module=$${dropout_structure_module} \
    --dropout_rates_filename=$${dropout_rates_filename} \
    --max_recycles=$${max_recycles} \
    --early_stop_tolerance=$${early_stop_tolerance} \
    --min_score=$${min_score} \
    --max_score=$${max_score} \
    --stop_recycling_below=$${stop_recycling_below} \
    --bfd_max_hits=$${bfd_max_hits} \
    --mgnify_max_hits=$${mgnify_max_hits} \
    --uniprot_max_hits=$${uniprot_max_hits} \
    --uniref_max_hits=$${uniref_max_hits} \
    --models_to_use=$${batch_model} \
    --start_prediction=$${batch_start} \
    --templates=$${templates} \
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --uniref30_database_path=$${uniref_database} \
    --pdb70_database_path=$${data_dir}/pdb70/pdb70

date
