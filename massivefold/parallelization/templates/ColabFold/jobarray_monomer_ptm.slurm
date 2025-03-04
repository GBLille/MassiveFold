
source $$(conda info --base)/etc/profile.d/conda.sh
conda activate mf-colabfold-1.5.5
#conda activate mf-colabfold

sequence_name=$sequence_name
run_name=$run_name
fafile=${input_dir}/$${sequence_name}.fasta
data_dir=$data_dir
scripts_dir=$scripts_dir

num_recycle=$num_recycle
recycle_early_stop_tolerance=$recycle_early_stop_tolerance
stop_at_score=$stop_at_score
BOOL_disable_cluster_profile=$disable_cluster_profile
BOOL_use_dropout=$use_dropout
if $$BOOL_disable_cluster_profile; then
  echo "Parameter --disable-cluster-profile set"
  disable_cluster_profile="--disable-cluster-profile"
fi
if $$BOOL_use_dropout; then
  echo "Parameter --use-dropout set"
  use_dropout="--use-dropout"
fi

output_dir=$$(realpath ${output_dir})
output_path="$${output_dir}/$${sequence_name}/$${run_name}/"
echo output_path is $$output_path
mkdir -vp $${output_path}/
path_to_msas=${output_dir}/$${sequence_name}/msas_colabfold

date

batches_file=${logs_dir}/${sequence_name}/${run_name}/$${sequence_name}_$${run_name}_batches.json
batch_model=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $$batches_file --element model)
batch_start=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $$batches_file --element start)
batch_end=$$($${scripts_dir}/get_batch.py --batch_id $$SLURM_ARRAY_TASK_ID --json_path $$batches_file --element end)

echo $$batch_model
num_models=$$(echo $$batch_model | cut -d "_" -f 2)
echo "using model $$num_models"
version_type="alphafold2_ptm"
model_type=$$version_type
echo "from version $$model_type"
num_seeds=$$(($$batch_end - $$batch_start + 1))
echo "predictions computed in the batch: $$num_seeds"

random_seed=$$(shuf -i 0-1000000 -n 1)
echo "Batch random seed: $$random_seed"

num_recycle=$num_recycle
recycle_early_stop_tolerance=$recycle_early_stop_tolerance
stop_at_score=$stop_at_score

echo "colabfold_batch
  $${path_to_msas}
  $${output_path}/batch_$${SLURM_ARRAY_TASK_ID}
  --data $${data_dir}
  --save-all
  --random-seed $$random_seed
  --num-seeds $$num_seeds
  --model-type $$model_type
  --model-order $$num_models
  --num-models 1
  --num-recycle $$num_recycle
  --recycle-early-stop-tolerance $$recycle_early_stop_tolerance
  --stop-at-score $$stop_at_score
  $$use_dropout
  $$disable_cluster_profile"

time colabfold_batch \
  $${path_to_msas} \
  $${output_path}/batch_$${SLURM_ARRAY_TASK_ID} \
  --data $${data_dir} \
  --save-all \
  --random-seed $$random_seed \
  --num-seeds $$num_seeds \
  --model-type $$model_type \
  --model-order $$num_models \
  --num-models 1 \
  --num-recycle $$num_recycle \
  --recycle-early-stop-tolerance $$recycle_early_stop_tolerance \
  --stop-at-score $$stop_at_score \
  $$use_dropout \
  $$disable_cluster_profile

date
