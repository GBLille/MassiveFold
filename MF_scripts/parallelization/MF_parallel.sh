#!/bin/bash

JOBNAME=test_multimer
prediction_number_per_model=3
batch_size=2
models_to_use=
cluster="jeanzay"

module load massivefold/1.0.0

./group_templates.py --cluster_name $cluster

# split the predictions in batches and store in json
./batching.py \
  --predictions_per_model=${prediction_number_per_model} \
  --batch_size=${batch_size} \
  --models_to_use=${models_to_use}

# Alignment job
./create_jobfile.py \
  --job_type alignment \
  --jobname $JOBNAME \
  --cluster_name $cluster

ALIGNMENT_ID=$(sbatch --parsable ${JOBNAME}_alignment.slurm)

# Inference by jobarray, waiting for alignment to end
./create_jobfile.py \
  --job_type jobarray \
  --jobname $JOBNAME \
  --cluster_name $cluster

ARRAY_ID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_ID ${JOBNAME}_jobarray.slurm)

#  Output organization and plots, waiting for inference to end
./create_jobfile.py \
  --job_type post_treatment \
  --jobname $JOBNAME \
  --cluster_name $cluster

sbatch --dependency=afterok:$ARRAY_ID ${JOBNAME}_post_treatment.slurm
