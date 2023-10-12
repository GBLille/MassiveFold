#!/bin/bash

#SBATCH --partition=base
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=200G
#SBATCH --gpus-per-node=1
#SBATCH --gres=gpu:1
#SBATCH -o ../log/%x.o%j
#SBATCH --array=0-$substitute_batch_number

date

module purge
module load massivefold/1.0.0

#for memory
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

#variables to adapt
data_dir=/shared/databases/2023-07-31
model_preset=multimer
max_template_date=2023-07-31
db_preset=full_dbs
use_gpu_relax=true
use_precomputed_msas=true
models_to_relax=best
alignments_only=false
dropout=false
dropout_rates_filename=
max_recycles=20
early_stop_tolerance=0.5
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
no_templates=false

jobname=$jobname
name=$$(basename -s .slurm $$jobname)
fafile=../input/$${name}.fasta

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
ln -s /home/"$$USER"/massivefold_runs/output_array/$${name}/msas $${output_dir}/$${name}

date

echo "massivefold-1.0.0
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

time massivefold-1.0.0 \
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
