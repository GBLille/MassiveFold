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
while test "$CONDA_PREFIX"
do
    conda deactivate
done

conda activate multicom-2.0.2

export PYTHONPATH="/home/Utils/MULTICOM3"

# for predictions on proteins that would typically be too long to fit into GPU memory
export TF_FORCE_UNIFIED_MEMORY="1"
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"

#fafile=../input/monomer.fasta
#flagfile=./flag/test_multimer.flg
#outputfile=../output

date

echo "time python /home/Utils/MULTICOM3/bin/monomer.py \
    --option_file=/home/Utils/MULTICOM3/bin/db_option \
    --fasta_path=../input/monomer.fasta \
    --run_img=True \
    --output_dir=../output/monomer_IMG"

time python /home/Utils/MULTICOM3/bin/monomer.py \
    --option_file=/home/Utils/MULTICOM3/bin/db_option \
    --fasta_path=../input/monomer.fasta \
    --run_img=True \
    --output_dir=../output/monomer_IMG

date


