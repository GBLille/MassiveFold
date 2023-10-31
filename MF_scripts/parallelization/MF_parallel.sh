#!/bin/bash

USAGE="\
./MF_parallel.sh -s arg -r arg -p arg -f arg [-b arg] [-m arg] [-c]\n\
./MF_parallel.sh -h for more details "

# help message
if [[ " ${@} " == *" -h "* ]] || [[ " ${@} " == *" --help "* ]]; then
  echo -e "\
Usage: $USAGE\n\
  Required arguments:\n\
    -s| --sequence: name of the sequence file without '.fasta'\n\
    -r| --run : name chosen for the run to store the outputs\n\
    -p| --predictions_per_model: number of predictions computed for each neural network model\n\
    -f| --parameters: path to the json file that contains the run's parameters\n\
\n\
  Facultative arguments:\n\
    -b| --batch_size: number of predictions per batch, default: 25\n\
    -m| --msas_precomputed: path to an output directory with msas already computed for the sequence
    -n| --select_model: Not available yet. Start a preliminary run, then select the top n models based on their ranking confidence median and launch another job exclusively with them.
\n\
  Facultative options:\n\
    -c| --calibrate: Does not work yet. Run a preleminary job to calibrate the batches size depending on the highest time for one prediction inference"
  exit 1
fi

calibration=false
predictions_per_model=67
batch_size=25

# argument parser
while true; do
  case "$1" in
    -s|--sequence)
      sequence_name=$2
      shift 2
      ;;
    -r|--run)
      run_name=$2
      shift 2
      ;;
    -p|--predictions_per_model)
      predictions_per_model=$2
      shift 2
      ;;
    -f|--parameters)
      parameters_file=$2
      shift 2
      ;;
    -b|--batch_size)
      batch_size=$2
      shift 2
      ;;
    -c|--calibrate)
      calibration=true
      shift
      ;;
    -m|--msas_precomputed)
      msas_precomputed=$2
      shift 2
      ;;
    -w|--wait_for_id)
      WAIT_FOR_ID=$2
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

# check mandatory args
if
  [ -z "$sequence_name" ] ||
  [ -z "$run_name" ] ||
  [ -z "$parameters_file" ]; then
  echo -e "Usage: $USAGE"
  exit 1
fi

echo "Run $run_name on sequence $sequence_name with $predictions_per_model predictions per model"

# launch calibration job and run post calibration MF_parallel waiting for it to end
if $calibration ; then
  echo -e "Running a priliminary job to calibrate the batches size.\n"
  echo "Feature is not ready, exiting"
  exit 1

  ./group_templates.py --parameters $parameters_file
  ./batching.py \
    --batch_size=5 \
    --predictions_per_model=5 \
    --sequence_name=$sequence_name \
    --run_name=calibration

  ./create_jobfile.py \
    --job_type=all \
    --sequence_name=${sequence_name} \
    --run_name=calibration \
    --path_to_parameters=${parameters_file}

  ALIGNMENT_ID=$(sbatch --parsable ${sequence_name}_calibration_alignment.slurm)
  ARRAY_ID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_ID ${sequence_name}_calibration_jobarray.slurm)
  END_CALIBRATION=$(sbatch --dependency=afterok:$ARRAY_ID ${sequence_name}_calibration_post_treatment.slurm)
  mkdir -p ../log_parallel/${sequence_name}/calibration/
  mv ${sequence_name}_calibration_* ../log_parallel/${sequence_name}/calibration/

  echo "Calibration not available yet, exiting."
  ./MF_parallel.sh \
    --sequence $sequence_name \
    --run $run_name \
    --predictions_per_model $predictions_per_model \
    --batch_size $batch_size \
    --parameters $parameters_file \
    --wait_for_id $END_CALIBRATION \
else
  echo "No calibration for the batch size."
fi


# Massivefold

module load massivefold/1.0.0

./group_templates.py --parameters $parameters_file

# split the predictions in batches and store in json
./batching.py \
  --predictions_per_model=${predictions_per_model} \
  --batch_size=${batch_size} \
  --models_to_use=${models_to_use} \
  --sequence_name=${sequence_name} \
  --run_name=${run_name}


if [ -z ${msas_precomputed} ]; then
  # Create and start alignment job
  ./create_jobfile.py \
    --job_type=alignment \
    --sequence_name=${sequence_name} \
    --run_name=${run_name} \
    --path_to_parameters=${parameters_file}

  ALIGNMENT_ID=$(sbatch --parsable ${sequence_name}_${run_name}_alignment.slurm)
elif [ -d  $msas_precomputed/msas ]; then
  echo "Using precomputed msas at $msas_precomputed"
  mkdir -p ../output_array/${sequence_name}/
  ln -s $(realpath $msas_precomputed/msas) ../output_array/${sequence_name}/
else
  echo "Directory $msas_precomputed does not exit or does not contain msas."
  exit 1
fi

# Create and launch inference jobarray 
./create_jobfile.py \
  --job_type=jobarray \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file} 

# Wait for alignment or calibration job, if neither needed, run
if ! [ -d $msas_precomputed/msas ]; then
  ARRAY_ID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_ID ${sequence_name}_${run_name}_jobarray.slurm)
elif 
  [ ! -z $WAIT_FOR_ID ] &&
  [ $(squeue -j $WAIT_FOR_ID | wc -l) -gt 1 ]; then
  ARRAY_ID=$(sbatch --parsable --dependency=afterok:$WAIT_FOR_ID ${sequence_name}_${run_name}_jobarray.slurm)
elif [! -z $WAIT_FOR_ID ]; then
  echo "Job $WAIT_FOR_ID didn't run properly or does not exist, exiting."
  exit 1
else
  ARRAY_ID=$(sbatch --parsable ${sequence_name}_${run_name}_jobarray.slurm)
fi

#  Create and start post treatment (output organization and plots)
./create_jobfile.py \
  --job_type=post_treatment \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file}

# Waiting for inference to end
sbatch --dependency=afterok:$ARRAY_ID ${sequence_name}_${run_name}_post_treatment.slurm

# Store jobiles and batches elements in logs

mkdir -p ../log_parallel/${sequence_name}/${run_name}/
mv ${sequence_name}_${run_name}_* ../log_parallel/${sequence_name}/${run_name}/

