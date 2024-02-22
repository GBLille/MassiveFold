date
# for memory
export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

#variables to adapt
sequence_name=$sequence_name
run_name=$run_name
echo "&${run_massivefold}&"
fafile=${input_dir}/$${sequence_name}.fasta
run_massivefold=$run_massivefold
data_dir=$data_dir
scripts_dir=$scripts_dir

output_dir=$$(realpath ${output_dir})
output_path="$${output_dir}/$${sequence_name}/$${run_name}/batch_$${SLURM_ARRAY_TASK_ID}"
echo output_path is $$output_path
mkdir -vp $${output_path}/$${sequence_name}
ln -s $${output_dir}/$${sequence_name}/msas $${output_path}/$${sequence_name}

# fixed variables
max_template_date=2023-07-31
use_precomputed_msas=true
alignments_only=false


# variable as MF_parallel parameters
model_preset=$AFM_run_model_preset
db_preset=$AFM_run_db_preset
use_gpu_relax=$AFM_run_use_gpu_relax
models_to_relax=$AFM_run_models_to_relax
dropout=$AFM_run_dropout
dropout_structure_module=$AFM_run_dropout_structure_module
dropout_rates_filename=$AFM_run_dropout_rates_filename
max_recycles=$AFM_run_max_recycles
early_stop_tolerance=$AFM_run_early_stop_tolerance
bfd_max_hits=$AFM_run_bfd_max_hits
mgnify_max_hits=$AFM_run_mgnify_max_hits
uniprot_max_hits=$AFM_run_uniprot_max_hits
uniref_max_hits=$AFM_run_uniref_max_hits
templates=$AFM_run_templates
min_score=$AFM_run_min_score
max_score=$AFM_run_max_score
stop_recycling_below=$AFM_run_stop_recycling_below

# parameters for each task in the job array
batches_path=${logs_dir}/${sequence_name}/${run_name}
batch_model=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element model)
batch_start=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element start)
batch_end=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}/$${sequence_name}_$${run_name}_batches.json --element end)

echo Model used is $$batch_model, batch computed between prediction $$batch_start and $$batch_end

echo "$${run_massivefold}
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
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta
"

time $${run_massivefold} \
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
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta
    
date
