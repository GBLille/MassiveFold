#!/bin/bash

# slurm parameters depending on your config 

date

# for memory
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

# fixed variables
run_alphafold_path= #e.g.:/gpfslocalsup/pub/anaconda-py3/2023.03/envs/MassiveFold-1.0.0/bin/run_alphafold.py on jean zay cluster
data_dir= #e.g.:/shared/databases/2023-07-31
max_template_date=2023-07-31
use_precomputed_msas=true
alignments_only=false

# variables to adapt
current_dir=$$(realpath .)
jobname=$jobname
name=$$(basename -s .slurm $$jobname)
fafile=../input/$${name}.fasta
model_preset=$substitute_model_preset
db_preset=$substitute_db_preset
use_gpu_relax=$substitute_use_gpu_relax
models_to_relax=$substitute_models_to_relax
dropout=$substitute_dropout
dropout_rates_filename=$substitute_dropout_rates_filename
max_recycles=$substitute_max_recycles
early_stop_tolerance=$substitute_early_stop_tolerance
bfd_max_hits=$substitute_bfd_max_hits
mgnify_max_hits=$substitute_mgnify_max_hits
uniprot_max_hits=$substitute_uniprot_max_hits
uniref_max_hits=$substitute_uniref_max_hits
no_templates=$substitute_no_templates

# parameters for each task in the job array
models_json=batches.json
batch_model=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element model)
batch_start=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element start)
batch_end=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element end)

models_to_use=$$batch_model

echo Model used is $$batch_model, batch computed between prediction $$batch_start and $$batch_end

output_dir="../output_array/$${name}/batch_$${SLURM_ARRAY_TASK_ID}_$${name}"
echo output_dir is $$output_dir
mkdir -vp $${output_dir}/$${name}
echo
ln -s $$current_dir/../output_array/$${name}/msas $${output_dir}/$${name}

date

echo "python3 $${run_alphafold_path}
    --fasta_paths=$${fafile}
    --output_dir=$${output_dir}
    --data_dir=$${data_dir}
    --db_preset=$${db_preset}
    --model_preset=$${model_preset}
    --models_to_relax=$${models_to_relax}
    --use_gpu_relax=$${use_gpu_relax}
    --max_template_date=$${max_template_date}
    --use_precomputed_msas=$${use_precomputed_msas}
    --num_predictions_per_model=$${batch_end}
    --alignments_only=$${alignments_only}
    --dropout=$${dropout}
    --dropout_rates_filename=$${dropout_rates_filename}
    --max_recycles=$${max_recycles}
    --early_stop_tolerance=$${early_stop_tolerance}
    --bfd_max_hits=$${bfd_max_hits}
    --mgnify_max_hits=$${mgnify_max_hits}
    --uniprot_max_hits=$${uniprot_max_hits}
    --uniref_max_hits=$${uniref_max_hits}
    --models_to_use=$${models_to_use}
    --start_prediction=$${batch_start}
    --no_templates=$${no_templates}
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta"

time python3 $${run_alphafold_path} \
    --fasta_paths=$${fafile} \
    --output_dir=$${output_dir} \
    --data_dir=$${data_dir} \
    --db_preset=$${db_preset} \
    --model_preset=$${model_preset} \
    --models_to_relax=$${models_to_relax} \
    --use_gpu_relax=$${use_gpu_relax} \
    --max_template_date=$${max_template_date} \
    --use_precomputed_msas=$${use_precomputed_msas} \
    --num_predictions_per_model=$${batch_end} \
    --alignments_only=$${alignments_only} \
    --dropout=$${dropout} \
    --dropout_rates_filename=$${dropout_rates_filename} \
    --max_recycles=$${max_recycles} \
    --early_stop_tolerance=$${early_stop_tolerance} \
    --bfd_max_hits=$${bfd_max_hits} \
    --mgnify_max_hits=$${mgnify_max_hits} \
    --uniprot_max_hits=$${uniprot_max_hits} \
    --uniref_max_hits=$${uniref_max_hits} \
    --models_to_use=$${models_to_use} \
    --start_prediction=$${batch_start} \
    --no_templates=$${no_templates} \
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta

date
