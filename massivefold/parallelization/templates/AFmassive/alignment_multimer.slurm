source $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-afmassive-1.1.5

# variables to adapt: used for alignment
sequence_name=$sequence_name
run_name=$run_name
run_massivefold=$run_massivefold
output=$output_dir
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir
uniref_database=$uniref_database

if [ -z "$$uniref_database" ]; then
  echo "No specified path for uniref30, use default ${data_dir}/uniref30/UniRef30_2021_03"
  uniref_database=$${data_dir}/uniref30/UniRef30_2021_03
fi

max_template_date=$max_template_date
#max_template_date=2023-07-31

model_preset=multimer
alignments_only=true
db_preset=full_dbs
bfd_max_hits=100000
mgnify_max_hits=501
uniprot_max_hits=50000
uniref_max_hits=10000
use_gpu_relax=false

echo "$${run_massivefold}
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
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat 
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt 
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt 
    --uniref30_database_path=$${uniref_database}
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta"


time $${run_massivefold} \
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
    --template_mmcif_dir=$${data_dir}/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=$${data_dir}/pdb_mmcif/obsolete.dat \
    --bfd_database_path=$${data_dir}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --pdb_seqres_database_path=$${data_dir}/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path=$${uniref_database} \
    --uniprot_database_path=$${data_dir}/uniprot/uniprot.fasta

date
