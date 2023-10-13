#!/usr/bin/env bash
##SBATCH --job-name=MF_inference
#SBATCH --account=$jeanzay_account

#SBATCH --error=../log_array/%x.%j
#SBATCH --output=../log_array/%x.%j

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --gpus-per-node=1
#SBATCH --array=0-$substitute_batch_number
#SBATCH --time=19:50:00
##SBATCH --qos=qos_gpu-t4         # Uncomment for job requiring more than 20h (max 16 GPUs)
#SBATCH -C $jeanzay_full_gpu             # Use gpu

module purge
#module load cpuarch/amd
module load massivefold/1.0.0

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH
export ALPHAFOLDDB=/gpfsdswork/dataset/AlphaFold-2.3.1

#variables to adapt
jobname=$jobname
name=$$(basename -s .slurm $$jobname)
fafile=../input/$${name}.fasta


output_dir="../output_array/$${name}/batch_$${SLURM_ARRAY_TASK_ID}_$${name}"
echo output_dir is $$output_dir
mkdir -vp $${output_dir}/$${name}

current_dir=$$(realpath .)
ln -s $$current_dir/../output_array/$${name}/msas $${output_dir}/$${name}

data_dir=$${ALPHAFOLDDB}/model_parameters/2.3.1
model_preset=multimer
max_template_date=2023-07-31
db_preset=full_dbs
use_gpu_relax=true
use_precomputed_msas=true
models_to_relax=best
alignments_only=false
dropout=true
dropout_rates_filename=
max_recycles=21
early_stop_tolerance=0.5
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
no_templates=true

models_json=batches.json
batch_model=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element model)
batch_start=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element start)
batch_end=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --element end)

models_to_use=$$batch_model
echo Model used is $$batch_model, batch computed between prediction $$batch_start and $$batch_end

date

echo "python3 $$(which run_alphafold.py)
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
    --uniref90_database_path=$${ALPHAFOLDDB}/uniref90/uniref90.fasta
    --mgnify_database_path=$${ALPHAFOLDDB}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${ALPHAFOLDDB}/pdb_mmcif
    --obsolete_pdbs_path=$${ALPHAFOLDDB}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${ALPHAFOLDDB}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb_seqres_database_path=$${ALPHAFOLDDB}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${ALPHAFOLDDB}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${ALPHAFOLDDB}/uniprot/uniprot.fasta"

time python3 $$(which run_alphafold.py) \
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
    --uniref90_database_path=$${ALPHAFOLDDB}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${ALPHAFOLDDB}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${ALPHAFOLDDB}/pdb_mmcif \
    --obsolete_pdbs_path=$${ALPHAFOLDDB}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${ALPHAFOLDDB}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${ALPHAFOLDDB}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${ALPHAFOLDDB}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${ALPHAFOLDDB}/uniprot/uniprot.fasta

date

