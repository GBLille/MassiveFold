#!/bin/bash

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