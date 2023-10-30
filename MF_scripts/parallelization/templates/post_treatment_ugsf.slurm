#!/bin/bash

#SBATCH --partition base
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=40G
##SBATCH --gpus-per-node=1
##SBATCH --gres=gpu:1
#SBATCH -o ../log/%x.o%j

date

source /ugsf/software/miniconda3/etc/profile.d/conda.sh
while test "$$CONDA_PREFIX"
do
    conda deactivate
done

conda activate massivefold-1.0.0
module load massivefold/1.0.0
date

sequence_name=$jobname
output=../output_array

# rename and move files
echo "time ./organize_outputs.py --batches_path $$output/$$sequence_name"
time ./organize_outputs.py --batches_path $$output/$$sequence_name

# delete not needed files
rm -r $$output/$$sequence_name/batch_*

# plots
echo "Generating plots"
echo "massivefold-1.0.0_plots --input_path $${output}/$${sequence_name} --plot_type for_each --top_n_predictions 5 --chosen_plots coverage,DM_plddt_PAE"

time massivefold-1.0.0_plots --input_path $${output}/$${sequence_name} --plot_type for_each --top_n_predictions 5  --chosen_plots coverage,DM_plddt_PAE

date
