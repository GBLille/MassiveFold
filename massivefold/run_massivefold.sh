#!/bin/bash


USAGE="\
./run_massivefold.sh -s str -r str -p int -f str [ -b int | [[-C str | -c] [-w int]] ] [-m str] [-n str] [-a] [-o]\n\
./run_massivefold.sh -h for more details "

# help message
if [[ " ${@} " == *" -h "* ]] || [[ " ${@} " == *" --help "* ]]; then
  echo -e "\
Usage: $USAGE\n\
  Required arguments:\n\
    -s| --sequence: path of the sequence(s) to infer, should be a 'fasta' file \n\
    -r| --run: name chosen for the run to organize in outputs.\n\
    -p| --predictions_per_model: number of predictions computed for each neural network model.\n\
    -f| --parameters: json file's path containing the parameters used for this run.\n\
\n\
  Facultative arguments:\n\
    -b| --batch_size: number of predictions per batch, should not be higher than -p (default: 25).\n\
    -m| --msas_precomputed: path to directory that contains computed msas.\n\
    -n| --top_n_models: uses the n neural network models with best ranking confidence from this run's path.\n\
    -w| --wall_time: total time available for calibration computations, unit is hours (default: 20).\n\
    -C| --calibration_from: path of a previous run to calibrate the batch size from (see --calibrate).\n\
\n\
  Facultative options:\n\
    -c| --calibrate: calibrate --batch_size value. Searches from the previous runs for the same 'fasta' path given\n\
        in --sequence and uses the longest prediction time found to compute the maximal number of predictions per batch.\n\
        This maximal number depends on the total time given by --wall_time.\n\
    -a| --recompute_msas: purges previous alignment step and recomputes msas.\n\
    -o| --only_msas: only compute alignments, the first step of MassiveFold."

  exit 1
fi

#module load massivefold/1.0.0
#conda activate massivefold-1.0.0

# default params
calibration=false
predictions_per_model=67
batch_size=25
wall_time=20
force_msas_computation=false
only_msas=false

# argument parser
while true; do
  case "$1" in
    -s|--sequence)
      sequence_file=$2
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
      path_to_run=$2
      shift 2
      ;;
    -w|--wall_time)
      wall_time=$2
      shift 2
      ;;
    -o|--only_msas)
      only_msas=true
      force_msas_computation=true
      shift
      ;;
    -a|--recompute_msas)
      force_msas_computation=true
      shift
      ;;
    *)
      break
      ;;
  esac
done

# check mandatory args
if
  [ -z "$sequence_file" ] ||
  [ -z "$run_name" ] ||
  [ -z "$parameters_file" ]; then
  echo -e "Usage: $USAGE"
  exit 1
fi

output_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['output_dir'])")
logs_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['logs_dir'])")
scripts_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['scripts_dir'])")

sequence_name=$(basename -s .fasta $sequence_file)

if [ ! -f ${sequence_file} ]; then
  echo "No sequence named ${sequence_name}.fasta in input directory $(dirname $sequence_file), exiting."
  exit 1
fi

# avoid overwriting run with a same name
# add an indicator (iteration number) if the run already exists 
if [ -d ${output_dir}/$sequence_name/$run_name ]; then
  echo "Run $run_name for $sequence_name already exists at ${output_dir}/$sequence_name/$run_name."
  echo "Starting new iteration of this run."
  i=1
  iteration=$run_name
  while [ -d ${output_dir}/$sequence_name/$iteration ]; do
    let i++
    iteration=${run_name}_$i
    echo "Trying $iteration"
  done
  run_name=$iteration
  echo -e "Current run is ${run_name}.\n"
fi

# calibration: check time taken for a previous run

number_of_runs=$(ls -l ${logs_dir}/${sequence_name} | wc -l)
if ! $calibration && [ -z $calibration_path ]; then
  echo "No calibration for the batch size."
elif $calibration && [ ! -z $calibration_path ]; then
  echo -e "Use either -c or -C, not both, exiting."
  exit 1
elif [ ! -z $calibration_path ] && [ ! -d $calibration_path ]; then
  echo "$calibration_path does not exist, exiting."
  exit 1
elif $calibration; then
  if [ ! -d ${logs_dir}/${sequence_name} ] && [ $number_of_runs -eq 0 ]; then
    echo "${sequence_name} has never been run, do try without calibration option, exiting."
    exit 1
  fi 
  echo "Calibrating this run's batch size."
  all_runs=$(find "${logs_dir}/$sequence_name" -mindepth 1 -maxdepth 1 -type d)
  echo "Searching for a completed preliminary run for $sequence_name"
  # search for completed run
  for run in $all_runs;
  do
    pred_nb=$(
    ${scripts_dir}/examine_run.py \
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
  ${scripts_dir}/examine_run.py \
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
    ${scripts_dir}/examine_run.py \
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
${scripts_dir}/batching.py \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --predictions_per_model=${predictions_per_model} \
  --batch_size=${batch_size} \
  --models_to_use=${models_to_use} \
  --path_to_parameters=${parameters_file}

# in case jobarrays start before the end of the script
mkdir -p ${logs_dir}/${sequence_name}/${run_name}/
cp ${sequence_name}_${run_name}_batches.json ${logs_dir}/${sequence_name}/${run_name}/


# align when forcing or no precomputed and detected msas

if [ ! -z $msas_precomputed ]; then
  echo "Using precomputed msas at $msas_precomputed"
elif [ -d ${output_dir}/${sequence_name}/msas/ ]; then
  echo -e "Detected msas for ${sequence_name} at ${output_dir}/${sequence_name}/msas/, \
  using them.\n"
  msas_precomputed="${output_dir}/${sequence_name}"
fi

waiting_for_alignment=false
conditions_to_align="[[ \$force_msas_computation = true ]] || \
                     ( [[ ! -d \${output_dir}/\${sequence_name}/msas/ ]] && \
                       [[ -z \$msas_precomputed ]] )"

if eval $conditions_to_align; then
  echo "Running alignment for $sequence_name"
  ${scripts_dir}/create_jobfile.py \
  --job_type=alignment \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file}

  ALIGNMENT_ID=$(sbatch --parsable ${sequence_name}_${run_name}_alignment.slurm)
  waiting_for_alignment=true
  if $only_msas; then
    echo "Only run sequence alignment."
    exit 1
  fi
elif [ -d  $msas_precomputed/msas ]; then
  echo "$msas_precomputed are valid."
  mkdir -p ${output_dir}/${sequence_name}/
  ln -s $(realpath $msas_precomputed/msas) ${output_dir}/${sequence_name}/
else
  echo "Directory $msas_precomputed does not exit or does not contain msas."
  exit 1
fi

# Create and launch inference jobarray 
${scripts_dir}/create_jobfile.py \
  --job_type=jobarray \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file} 

# Only wait for alignment if not precomputed
if [[ $waiting_for_alignment = true  ]]; then
  ARRAY_ID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_ID ${sequence_name}_${run_name}_jobarray.slurm)
else
  ARRAY_ID=$(sbatch --parsable ${sequence_name}_${run_name}_jobarray.slurm)
fi

# Create and start post treatment (output organization axnd plots)
# Waiting for inference to end
${scripts_dir}/create_jobfile.py \
  --job_type=post_treatment \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file}

sbatch --dependency=afterok:$ARRAY_ID ${sequence_name}_${run_name}_post_treatment.slurm

# Store jobiles and batches elements in logs
mkdir -p ${logs_dir}/${sequence_name}/${run_name}/
mv ${sequence_name}_${run_name}_* ${logs_dir}/${sequence_name}/${run_name}/
