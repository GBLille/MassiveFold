#!/usr/bin/env bash

#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks per node
#SBATCH --cpus-per-task=10    # Number of OpenMP threads per task
##SBATCH --gpus-per-node=1    # Number of GPUs per node
##SBATCH -C a100              # Use A100 partition
#SBATCH --hint=nomultithread # Disable hyperthreading
##SBATCH --job-name=MF_alignment # Jobname
#SBATCH --output=../log/%x.o%j      # Output file %x is the jobname, %j the jobid
#SBATCH --error=../log/%x.o%j       # Error file
#SBATCH --time=10:00:00      # Expected runtime HH:MM:SS (max 100h)
#SBATCH --account=fvp@cpu #
##
## Please, refer to comments below for
## more information about these 4 last options.
##SBATCH --account=<account>@a100      # To specify cpu accounting: <account> = echo $$IDRPROJ
##SBATCH --partition=<partition>       # To specify partition (see IDRIS web site for more info)
##SBATCH --qos=qos_gpu-dev             # Uncomment for job requiring less than 2 hours
##SBATCH --qos=qos_gpu-t4              # Uncomment for job requiring more than 20h (max 16 GPUs)

module purge
module load cpuarch/amd
#module load alphafold/2.3.1
module purge
module load massivefold/1.0.0

export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH
export ALPHAFOLDDB=/gpfsdswork/dataset/AlphaFold-2.3.1


jobname=$jobname
fafile=../input/$${jobname}.fasta
output=../output_array/

data_dir=$${ALPHAFOLDDB}/model_parameters/2.3.1
database_dir=$${ALPHAFOLDDB}
model_preset=multimer
max_template_date=2023-07-31
db_preset=full_dbs
num_predictions_per_model=5
alignments_only=true
dropout=false
dropout_rates_filename=
max_recycles=20
early_stop_tolerance=0.5
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
models_to_use=
start_prediction=1
no_templates=false
models_to_relax=none
use_gpu_relax=false
use_precomputed_msas=false

date

echo "python3 $$(which run_alphafold.py)
    --fasta_paths=$${fafile}
    --output_dir=$${output}
    --data_dir=$${data_dir}
    --db_preset=$${db_preset}
    --model_preset=$${model_preset}
    --models_to_relax=$${models_to_relax}
    --use_gpu_relax=$${use_gpu_relax}
    --max_template_date=$${max_template_date}
    --use_precomputed_msas=$${use_precomputed_msas}
    --num_predictions_per_model=$${num_predictions_per_model}
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
    --start_prediction=$${start_prediction}
    --no_templates=$${no_templates}
    --uniref90_database_path=$${database_dir}/uniref90/uniref90.fasta
    --mgnify_database_path=$${database_dir}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${database_dir}/pdb_mmcif
    --obsolete_pdbs_path=$${database_dir}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${database_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb_seqres_database_path=$${database_dir}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${database_dir}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${database_dir}/uniprot/uniprot.fasta"

time python3 $$(which run_alphafold.py) \
    --fasta_paths=$${fafile} \
    --output_dir=$${output} \
    --data_dir=$${data_dir} \
    --db_preset=$${db_preset} \
    --model_preset=$${model_preset} \
    --models_to_relax=$${models_to_relax} \
    --use_gpu_relax=$${use_gpu_relax} \
    --max_template_date=$${max_template_date} \
    --use_precomputed_msas=$${use_precomputed_msas} \
    --num_predictions_per_model=$${num_predictions_per_model} \
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
    --start_prediction=$${start_prediction} \
    --no_templates=$${no_templates} \
    --uniref90_database_path=$${database_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${database_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${database_dir}/pdb_mmcif \
    --obsolete_pdbs_path=$${database_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${database_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${database_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${database_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${database_dir}/uniprot/uniprot.fasta

date
