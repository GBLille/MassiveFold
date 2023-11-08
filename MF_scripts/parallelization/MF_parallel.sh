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
    -n|--select_model_from_run)
      path_to_run="$2"
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

if $calibration ; then
  echo -e "Running a preliminary job to calibrate the batches size.\n"
  echo "Calibration not available yet, exiting."
  exit 1
else
  echo "No calibration for the batch size."
fi

: '
if [ -n $path_to_run ]; then
  echo "Running with the 5 best models of the $path_to_run run on $sequence_name"
  exit 1
  echo "verify ../output_array/$sequence_name/$path_to_run/ranking_debug.json"
  if [ -f ../output_array/$sequence_name/$path_to_run/ranking_debug.json ]; then
    echo entered
    models_to_use=$(./analyze.py --log_dir ../log_parallel/$sequence_name/$path_to_run/ --get models)
    echo $models_to_use
    exit 1
  else
    echo "not entered"
    echo "Either $path_to_run isn't finished or it doesn't exist, exiting."
    exit 1
  fi
fi
'

# Massivefold

module load massivefold/1.0.0

./group_templates.py --parameters $parameters_file

# split the predictions in batches and store in json
./batching.py \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --predictions_per_model=${predictions_per_model} \
  --batch_size=${batch_size} \
  --models_to_use=${models_to_use} \
  --path_to_parameters=${parameters_file}

# in case jobarrays start before the end of the script
mkdir -p ../log_parallel/${sequence_name}/${run_name}/
cp ${sequence_name}_${run_name}_batches.json ../log_parallel/${sequence_name}/${run_name}/

# starts alignment only when not pre-existing
if [ -d ../output_array/${sequence_name}/msas/ ]; then
  echo "Detected msas for the same ${sequence_name} at ../output_array/${sequence_name}/msas/, using them."
elif [ -z ${msas_precomputed} ]; then
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

# Only wait for alignment if not precomputed
if ! [ -d $msas_precomputed/msas ]; then
  ARRAY_ID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_ID ${sequence_name}_${run_name}_jobarray.slurm)
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

