#!/bin/bash


USAGE="\
./MF_parallel.sh -s str -r str -p int -f str [-m str] [-n str] [[-b int] | [-C str] | [-c]]\n\
./MF_parallel.sh -h for more details "

# help message
if [[ " ${@} " == *" -h "* ]] || [[ " ${@} " == *" --help "* ]]; then
  echo -e "\
Usage: $USAGE\n\
  Required arguments:\n\
    -s| --sequence: name of the sequence file without '.fasta'.\n\
    -r| --run: name chosen for the run to store the outputs.\n\
    -p| --predictions_per_model: number of predictions computed for each neural network model.\n\
    -f| --parameters: path to the json file that contains the run's parameters.\n\
\n\
  Facultative arguments:\n\
    -b| --batch_size: number of predictions per batch, default: 25.\n\
    -m| --msas_precomputed: path to output folder containing already computed msas.\n\
    -n| --top_n_models: path of a completed run, use the 5 best models from the location.\n\
    -C| --calibration_from: path of a previous run to calibrate the batch size.\n\
    -
\n\
  Facultative options:\n\
    -c| --calibrate_batch_size: set the --batch_size by computing the maximal number of prediction per batch.
It searches for previous runs on the sequence and use the longest prediction time found."
  exit 1
fi

module load massivefold/1.0.0

# default params
calibration=false
predictions_per_model=67
batch_size=25
wall_time=20

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
    -m|--msas_precomputed)
      msas_precomputed=$2
      shift 2
      ;;
    -b|--batch_size)
      batch_size=$2
      shift 2
      ;;
    -c|--calibrate_batch_size)
      calibration=true
      shift
      ;;
    -C|--calibration_from)
      calibration_path=$2
      shift 2
      ;;
    -n|--top_n_model)
      path_to_run="$2"
      shift 2
      ;;
    -w|--wall_time)
      wall_time="$2"
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

# avoid overwriting run with a same name
# add an indicator (iteration number) if the run already exists 
if [ -d ../output_array/$sequence_name/$run_name ]; then
  echo "Run $run_name for $sequence_name already exists at ../output_array/$sequence_name/$run_name."
  echo "Starting new iteration of this run."
  i=1
  iteration=$run_name
  while [ -d ../output_array/$sequence_name/$iteration ]; do
    let i++
    iteration=${run_name}_$i
    echo "Trying $iteration"
  done
  run_name=$iteration
  echo -e "Current run is ${run_name}.\n"
fi

# check mandatory args
if
  [ -z "$sequence_name" ] ||
  [ -z "$run_name" ] ||
  [ -z "$parameters_file" ]; then
  echo -e "Usage: $USAGE"
  exit 1
fi


# calibration: check time taken for a previous run
number_of_runs=$(ls -l ../log_parallel/${sequence_name} | wc -l)
if ! $calibration && [ -z $calibration_path ]; then
  echo "No calibration for the batch size."
elif $calibration && [ ! -z $calibration_path ]; then
  echo -e "Use either -c or -C, not both, exiting."
  exit 1
elif [ ! -z $calibration_path ] && [ ! -d $calibration_path ]; then
  echo "$calibration_path does not exist, exiting."
  exit 1
elif $calibration; then
  if [ ! -d ../log_parallel/${sequence_name} ] && [ $number_of_runs -eq 0 ]; then
    echo "${sequence_name} has never been run, do try without calibration option, exiting."
    exit 1
  fi 
  echo "Calibrating this run's batch size."
  all_runs=$(find "../log_parallel/$sequence_name" -mindepth 1 -maxdepth 1 -type d)
  echo "Searching for a completed preliminary run for $sequence_name"
  # search for completed run
  for run in $all_runs;
  do
    pred_nb=$(
    ./examine_run.py \
      --get=batch_calibration \
      --input=${run} \
      --wall_time=${wall_time}
    )
    if [[ $pred_nb =~ ^[0-9]+$ ]]; then
      if [ -z $lowest_pred_nb ]; then
        lowest_pred_nb=$pred_nb
        echo "Found first batch size candidate: ${lowest_pred_nb} at ${run}"
      elif [ $pred_nb -le $lowest_pred_nb ]; then
        lowest_pred_nb=$pred_nb
        echo "Found new batch size candidate: ${lowest_pred_nb} at ${run}"
      fi
    fi
  done
  batch_size=$lowest_pred_nb 
  
  if [ -z $lowest_pred_nb ]; then
    echo "No preliminary run completed for $sequence_name, exiting."
    exit 1
  fi

elif [ ! -z $calibration_path ]; then
  if [ ! -d $calibration_path ]; then
    echo "$calibration_path does not exist or not completed, exiting."
    exit 1
  fi

  echo "Calibrating this run's batch size."
  pred_nb=$(
  ./examine_run.py \
    --get=batch_calibration \
    --input=${calibration_path} \
    --wall_time=${wall_time}
  )
  if [[ $pred_nb =~ ^[0-9]+$ ]]; then
    batch_size=$pred_nb
  fi
fi

if $calibration || [ ! -z $calibration_path ]; then
  echo "Number of prediction under wall time: ${batch_size}"
  if [ $batch_size -gt $predictions_per_model ]; then
    batch_size=$predictions_per_model
    echo "Adapting batch size according to -p ${predictions_per_model}"
  fi
  echo -e "Calibrated batch size: ${batch_size} \n"
fi


if [ ! -z $path_to_run ]; then
  echo "Running with the 5 best models of the run located at path_to_run."
  if [ -f ${path_to_run}/ranking_debug.json ]; then
    echo "Using ${path_to_run} run to evaluate the 5 best models."
    models_to_use=$(
    ./examine_run.py \
      --get=models \
      --input=${path_to_run}
    )
    echo -e "Using following models: $models_to_use \n"
  else
    echo "Either run $path_to_run is still running or does not exist, exiting."
    exit 1
  fi
fi

echo "Run $run_name on sequence $sequence_name with $predictions_per_model predictions per model"

# Massivefold

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
  echo -e "Detected msas for ${sequence_name} located ../output_array/${sequence_name}/msas/, using them.\n"
  msas_precomputed="../output_array/${sequence_name}"
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

#  Create and start post treatment (output organization axnd plots)
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
