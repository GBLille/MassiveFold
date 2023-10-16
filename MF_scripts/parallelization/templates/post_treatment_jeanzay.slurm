#!/bin/bash

#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of tasks per node
#SBATCH --cpus-per-task=10    # Number of OpenMP threads per task
##SBATCH --gpus-per-node=1    # Number of GPUs per node
##SBATCH -C $jeanzay_full_gpu              # Use A100 partition
#SBATCH --hint=nomultithread # Disable hyperthreading
##SBATCH --job-name=MF_post_treatment # Jobname
#SBATCH --output=../log/%x.%j.log
#SBATCH --error=../log/%x.%j.log
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

sequence_name=$jobname
output=../output_array

# rename and move files
echo "time ./organize_outputs.py --batches_path $$output/$$sequence_name"
time ./organize_outputs.py --batches_path $$output/$$sequence_name

# delete not needed files
rm -r $$output/$$sequence_name/batch_*

# plots
echo "Generating plots"

echo "python3 /gpfswork/rech/fvp/commun/plots/MF_plots.py --input_path $${output}/$${sequence_name} --plot_type for_each --top_n_predictions 5 --chosen_plots coverage,DM_plddt_PAE"
time python3 /gpfswork/rech/fvp/commun/plots/MF_plots.py --input_path $${output}/$${sequence_name} --plot_type for_each --top_n_predictions 5 --chosen_plots coverage,DM_plddt_PAE

echo "python3 /gpfswork/rech/fvp/commun/plots/MF_plots.py --input_path $${output}/$${sequence_name} --top_n_predictions 5 --chosen_plots CF_PAEs"
time python3 /gpfswork/rech/fvp/commun/plots/MF_plots.py --input_path $${output}/$${sequence_name} --top_n_predictions 5 --chosen_plots CF_PAEs

date
