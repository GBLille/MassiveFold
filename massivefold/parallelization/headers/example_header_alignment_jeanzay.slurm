#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --hint=nomultithread
#SBATCH --output=${logs_dir}/${sequence_name}/${run_name}/alignment_%j.log
#SBATCH --error=${logs_dir}/${sequence_name}/${run_name}/alignment_%j.log
#SBATCH --time=$jeanzay_alignment_time
#SBATCH --account=$jeanzay_project@cpu
#SBATCH --partition=prepost
##SBATCH --qos=qos_gpu-t4              # Uncomment for job requiring more than 20h (max 16 GPUs)

#tmp dir on Jean Zay
export TMP=$$JOBSCRATCH
export TMPDIR=$$JOBSCRATCH

