#!/bin/bash

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
##SBATCH --qos=qos_gpu-t4               # Uncomment for job requiring more than 20h (max 16 GPUs)
#SBATCH -C $jeanzay_gpu_with_memory     # Use gpu

#tmp dir on Jean Zay
export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH

# to run on H100:
#module purge
#module load arch/h100

# to run on A100:
#module purge
#module load arch/a100

#For memory, to activate for very large systems for AFmassive or ColabFold
#export TF_FORCE_UNIFIED_MEMORY='1'
#export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'
