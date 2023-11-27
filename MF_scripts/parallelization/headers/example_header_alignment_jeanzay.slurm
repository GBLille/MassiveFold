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
