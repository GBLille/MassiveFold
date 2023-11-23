#!/usr/bin/env bash
##SBATCH --job-name=MF_inference
#SBATCH --account=$jeanzay_account

#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/jobarray_%a.log

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --gpus-per-node=1
#SBATCH --array=0-$substitute_batch_number
#SBATCH --time=$jeanzay_jobarray_time
##SBATCH --qos=qos_gpu-dev             # Uncomment for job requiring less than 2 hours
##SBATCH --qos=qos_gpu-t4         # Uncomment for job requiring more than 20h (max 16 GPUs)
#SBATCH -C $jeanzay_gpu_with_memory             # Use gpu

module purge
module load massivefold/1.0.0

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH
export ALPHAFOLDDB=$data_dir

#variables to adapt
sequence_name=$sequence_name
run_name=$run_name
fafile=${input_dir}/$${sequence_name}.fasta
output_dir="${output_dir}/$${sequence_name}/$${run_name}/batch_$${SLURM_ARRAY_TASK_ID}"
echo output_dir is $$output_dir
mkdir -vp $${output_dir}/$${sequence_name}

use_precomputed_msas=true
alignments_only=false

#path_run_AF2=$$(which run_alphafold.py)
#path_run_AF2=/gpfs7kw/linkhome/rech/genasu01/urh88pc/git/MassiveFold_dev/run_alphafold.py
path_run_AF2=$run_massivefold


current_dir=$$(realpath .)
ln -s $$current_dir/${output_dir}/$${sequence_name}/msas $${output_dir}/$${sequence_name}
data_dir=$${ALPHAFOLDDB}/model_parameters/2.3.1

model_preset=$MF_run_model_preset
max_template_date=2023-07-31
db_preset=$MF_run_db_preset
use_gpu_relax=$MF_run_use_gpu_relax
models_to_relax=$MF_run_models_to_relax
dropout=$MF_run_dropout
dropout_structure_module=$MF_run_dropout_structure_module
dropout_rates_filename=$MF_run_dropout_rates_filename
max_recycles=$MF_run_max_recycles
early_stop_tolerance=$MF_run_early_stop_tolerance
bfd_max_hits=$MF_run_bfd_max_hits
mgnify_max_hits=$MF_run_mgnify_max_hits
uniprot_max_hits=$MF_run_uniprot_max_hits
uniref_max_hits=$MF_run_uniref_max_hits
templates=$MF_run_templates
score_threshold_output=$MF_run_score_threshold_output
max_score_stop=$MF_run_max_score_stop

batches_path=../log_parallel/${sequence_name}/${run_name}/
batch_model=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}./$${sequence_name}_$${run_name}_batches.json --element model)
batch_start=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}./$${sequence_name}_$${run_name}_batches.json --element start)
batch_end=$$(./get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $${batches_path}./$${sequence_name}_$${run_name}_batches.json --element end)

models_to_use=$$batch_model
echo Model used is $$batch_model, batch computed between prediction $$batch_start and $$batch_end

date

echo "python3 $${path_run_AF2}
    --fasta_paths=$${fafile}
    --output_dir=$${output_dir}
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
    --dropout_rates_filename=$${dropout_rates_filename}
    --dropout_structure_module=$${dropout_structure_module}
    --max_recycles=$${max_recycles}
    --early_stop_tolerance=$${early_stop_tolerance}
    --bfd_max_hits=$${bfd_max_hits}
    --mgnify_max_hits=$${mgnify_max_hits}
    --uniprot_max_hits=$${uniprot_max_hits}
    --uniref_max_hits=$${uniref_max_hits}
    --models_to_use=$${models_to_use}
    --start_prediction=$${batch_start}
    --templates=$${templates}
    --uniref90_database_path=$${ALPHAFOLDDB}/uniref90/uniref90.fasta
    --mgnify_database_path=$${ALPHAFOLDDB}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${ALPHAFOLDDB}/pdb_mmcif
    --obsolete_pdbs_path=$${ALPHAFOLDDB}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${ALPHAFOLDDB}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb_seqres_database_path=$${ALPHAFOLDDB}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${ALPHAFOLDDB}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${ALPHAFOLDDB}/uniprot/uniprot.fasta
    --score_threshold_output=$${score_threshold_output} 
    --max_score_stop=$${max_score_stop}
"

time python3 $${path_run_AF2} \
    --fasta_paths=$${fafile} \
    --output_dir=$${output_dir} \
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
    --dropout_rates_filename=$${dropout_rates_filename} \
    --dropout_structure_module=$${dropout_structure_module} \
    --max_recycles=$${max_recycles} \
    --early_stop_tolerance=$${early_stop_tolerance} \
    --bfd_max_hits=$${bfd_max_hits} \
    --mgnify_max_hits=$${mgnify_max_hits} \
    --uniprot_max_hits=$${uniprot_max_hits} \
    --uniref_max_hits=$${uniref_max_hits} \
    --models_to_use=$${models_to_use} \
    --start_prediction=$${batch_start} \
    --templates=$${templates} \
    --uniref90_database_path=$${ALPHAFOLDDB}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${ALPHAFOLDDB}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${ALPHAFOLDDB}/pdb_mmcif \
    --obsolete_pdbs_path=$${ALPHAFOLDDB}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${ALPHAFOLDDB}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${ALPHAFOLDDB}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${ALPHAFOLDDB}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${ALPHAFOLDDB}/uniprot/uniprot.fasta \
    --score_threshold_output=$${score_threshold_output} \
    --max_score_stop=$${max_score_stop}
    

date
