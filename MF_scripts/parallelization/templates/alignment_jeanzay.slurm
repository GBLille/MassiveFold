#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --output=../log_parallel/${sequence_name}/${run_name}/alignment_%j.log
#SBATCH --error=../log_parallel/${sequence_name}/${run_name}/alignment_%j.log
#SBATCH --time=$jeanzay_alignment_time
#SBATCH --account=$jeanzay_project@cpu
##SBATCH --qos=qos_gpu-t4              # Uncomment for job requiring more than 20h (max 16 GPUs)


module purge
module load massivefold/1.0.0

export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH
export ALPHAFOLDDB=/gpfsdswork/dataset/AlphaFold-2.3.1


run_name=$run_name
sequence_name=$sequence_name
fafile=../input/$${sequence_name}.fasta
output=../output_array/

run_AF_path=$$(which run_alphafold.py) 
run_AF_path=/gpfs7kw/linkhome/rech/genasu01/urh88pc/git/MassiveFold_dev/run_alphafold.py

data_dir=$${ALPHAFOLDDB}/model_parameters/2.3.1
database_dir=$${ALPHAFOLDDB}
model_preset=multimer
max_template_date=2023-07-31
db_preset=full_dbs
alignments_only=true
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
use_gpu_relax=false

date

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
    --uniref90_database_path=$${database_dir}/uniref90/uniref90.fasta 
    --mgnify_database_path=$${database_dir}/mgnify/mgy_clusters_2022_05.fa 
    --template_mmcif_dir=$${database_dir}/pdb_mmcif 
    --obsolete_pdbs_path=$${database_dir}/pdb_mmcif/obsolete.dat 
    --bfd_database_path=$${database_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt 
    --pdb_seqres_database_path=$${database_dir}/pdb_seqres/pdb_seqres.txt 
    --uniref30_database_path=$${database_dir}/uniref30/UniRef30_2021_03 
    --uniprot_database_path=$${database_dir}/uniprot/uniprot.fasta"


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
    --uniref90_database_path=$${database_dir}/uniref90/uniref90.fasta \
    --mgnify_database_path=$${database_dir}/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$${database_dir}/pdb_mmcif \
    --obsolete_pdbs_path=$${database_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${database_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${database_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${database_dir}/uniref30/UniRef30_2021_03 \
    --uniprot_database_path=$${database_dir}/uniprot/uniprot.fasta

date
