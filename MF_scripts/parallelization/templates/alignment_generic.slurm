#!/bin/bash

# slurm parameters to launch on cpu

# activate AF/MF conda env
# e.g: module load massivefold/1.0.0

# fixed variables
run_alphafold_path= #e.g.:'run_alphafold_path=/gpfslocalsup/pub/anaconda-py3/2023.03/envs/MassiveFold-1.0.0/bin/run_alphafold.py' for jean zay cluster
data_dir= #e.g.:'data_dir=/shared/databases/2023-07-31'
max_template_date=2023-07-31

output=../output_array/
use_precomputed_msas=false
alignments_only=true
max_template_date=2023-07-31

# variables to adapt: used for alignment
jobname=$jobname
fafile=../input/$$jobname.fasta
model_preset=$substitute_model_preset
db_preset=$substitute_db_preset
bfd_max_hits=$substitute_bfd_max_hits
mgnify_max_hits=$substitute_mgnify_max_hits
uniprot_max_hits=$ubstitute_uniprot_max_hits
uniref_max_hits=$substitute_uniref
no_templates=$substitute_no_templates

# variable to ignore: unused for alignment
num_predictions_per_model=5
dropout=false
dropout_rates_filename=
max_recycles=20
early_stop_tolerance=0.5
models_to_use=
start_prediction=1
models_to_relax=none
use_gpu_relax=false

date

echo "python3 $${run_alphafold_path}
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
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta"

time python3 ${run_alphafold_path} \
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
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta

date
