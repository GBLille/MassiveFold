#!/bin/bash

# slurm parameters to launch on cpu

# activate AF/MF conda env
# e.g: module load massivefold/1.0.0

# fixed variables
run_AF_path= #e.g.:'run_alphafold_path=/gpfslocalsup/pub/anaconda-py3/2023.03/envs/MassiveFold-1.0.0/bin/run_alphafold.py' for jean zay cluster
data_dir= #e.g.:'data_dir=/shared/databases/2023-07-31'
max_template_date=2023-07-31

# variables to adapt: used for alignment
run_name=$run_name
sequence_name=$sequence_name
fafile=../input/$${sequence_name}.fasta
output=../output_array/

model_preset=multimer
max_template_date=2023-07-31
db_preset=full_dbs
alignments_only=true
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
use_gpu_relax=false

echo "python3 $${run_AF_path} 
    --alignments_only=$${alignments_only} 
    --fasta_paths=$${fafile} 
    --output_dir=$${output} 
    --data_dir=$${data_dir} 
    --db_preset=$${db_preset} 
    --use_gpu_relax=$${use_gpu_relax} 
    --model_preset=$${model_preset} 
    --max_template_date=$${max_template_date} 
    --bfd_max_hits=$${bfd_max_hits} 
    --mgnify_max_hits=$${mgnify_max_hits} 
    --uniprot_max_hits=$${uniprot_max_hits} 
    --uniref_max_hits=$${uniref_max_hits} 
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta 
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa 
    --template_mmcif_dir=$${data_dir}/pdb_mmcif 
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat 
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt 
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt 
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03 
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta"


time python3 $${run_AF_path} \
    --alignments_only=$${alignments_only} \
    --fasta_paths=$${fafile} \
    --output_dir=$${output} \
    --data_dir=$${data_dir} \
    --db_preset=$${db_preset} \
    --use_gpu_relax=$${use_gpu_relax} \
    --model_preset=$${model_preset} \
    --max_template_date=$${max_template_date} \
    --bfd_max_hits=$${bfd_max_hits} \
    --mgnify_max_hits=$${mgnify_max_hits} \
    --uniprot_max_hits=$${uniprot_max_hits} \
    --uniref_max_hits=$${uniref_max_hits} \
    --uniref90_database_path=$${data_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${data_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${data_dir}/pdb_mmcif \
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${data_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta

date
