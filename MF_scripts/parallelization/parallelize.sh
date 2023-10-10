#!/bin/bash

JOBNAME=$1
prediction_number_per_model=
batch_size=
models_to_use=

#split the predictions in batches and store in json
./batching.py --predictions_per_model ${prediction_number_per_model} --models_to_use ${models_to_use} --batch_size ${batch_size}

#create the slurm job array
./create_jobfile.py --jobname $JOBNAME

#launch the job array and get its id
ARRAY_ID=$(sbatch --parsable ${JOBNAME}.slurm)

#Job: organize output and make plots
./create_jobfile.py --with_output --jobname $JOBNAME

#wait for previous job to end to launch the job
sbatch --dependency=afterok:$ARRAY_ID order_and_plots.slurm
