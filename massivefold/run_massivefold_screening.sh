#!/bin/bash

USAGE="\
./run_massivefold_screening.sh -s str -r str -p int -f str -t str [ -b int | [[-C str | -c] [-w int]] ] [-m str] [-n str] [-a] [-o]\n\
./run_massivefold_screening.sh -h for more details "

# help message
if [[ " ${@} " == *" -h "* ]] || [[ " ${@} " == *" --help "* ]]; then
  echo -e "\
Usage: $USAGE\n\
  Required arguments:\n\
    -s| --sequence: path of the fasta file containing sequence(s) used for screening.\n\
    -l| --ligands: csv file containing the list of ligands to use for screening. \n\
    -f| --parameters: json file's path containing the parameters used for the screening.\n\
\n\
  Facultative arguments:\n\
    -p| --predictions_per_model: number of seed used. Each seed will have 5 samples predicted.\n\
        In total, with -p n, you will have 5n predictions computed.\n\
    -m| --msas_precomputed: path to directory that contains computed msas.\n\
    -j| --jobid: jobid of an alignment job to wait for inference, skips the alignments.\n\
\n\
  Facultative options:\n\
    -o| --only_msas: only compute alignments, the first step of MassiveFold"
  exit 1
fi

# default params
calibration=false
predictions_per_model=1
batch_size=1
only_msas=false

# argument parser
while true; do
  case "$1" in
    -s|--sequence)
      sequence_file=$2
      shift 2
      ;;
    -l|--ligands)
      ligand_file=$2
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
    -o|--only_msas)
      only_msas=true
      shift
      ;;
    -j|--jobid)
      wait_for_jobid=$2
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

# check mandatory args
if
  [ -z "$sequence_file" ] ||
  [ -z "$ligand_file" ] ||
  [ -z "$parameters_file" ]; then
  echo -e "Usage: $USAGE"
  echo "Argument(s) missing"
  exit 1
fi

output_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['output_dir'])")
logs_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['logs_dir'])")
scripts_dir=$(cat $parameters_file | python3 -c "import sys, json; print(json.load(sys.stdin)['massivefold']['scripts_dir'])")

sequence_name=$(basename -s .fasta $sequence_file)
ligand_filename=$(basename -s .csv $ligand_file)
run_name=$ligand_filename

if [ ! -f ${sequence_file} ]; then
  echo "No sequence named ${sequence_name}.fasta in input directory $(dirname $sequence_file), exiting."
  exit 1
fi
if [ ! -f ${ligand_file} ]; then
  echo "No sequence named ${ligand_filename}.csv in input directory $(dirname $ligand_file), exiting."
  exit 1
fi

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

echo "Run $run_name on sequence $sequence_name with $predictions_per_model predictions per model"

# Massivefold

# split the predictions in batches and store in json
${scripts_dir}/batching.py \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --predictions_per_model=${predictions_per_model} \
  --batch_size=${batch_size} \
  --models_to_use=${models_to_use} \
  --path_to_parameters=${parameters_file} \
  --to_screen $ligand_file \
  --tool AlphaFold3 \
  || { echo "Dividing predictions into batches failed. Exiting."; exit 1; }

# in case jobarrays start before the end of the script
mkdir -p ${logs_dir}/${sequence_name}/${run_name}/
cp ${sequence_name}_${run_name}_batches.json ${logs_dir}/${sequence_name}/${run_name}/

# align when forcing or no precomputed and detected msas
if [ ! -z $msas_precomputed ]; then
  echo "Using precomputed msas at $msas_precomputed"
elif [ -d ${output_dir}/${sequence_name}/msas_alphafold3/ ]; then
  echo -e "Detected msas compatible with af3 for ${sequence_name} at ${output_dir}/${sequence_name}/msas_alphafold3/, using them.\n"
  msas_precomputed="${output_dir}/${sequence_name}/msas_alphafold3"
fi

waiting_for_alignment=false

conditions_to_align="( [[ ! -d \${output_dir}/\${sequence_name}/msas_alphafold3/ ]] && [[ -z \$msas_precomputed ]] )"

using_jobid=false
if [ ! -z $wait_for_jobid ]; then
  echo "Waiting for alignment job $wait_for_jobid"
  ALIGNMENT_ID=$wait_for_jobid
  using_jobid=true
  waiting_for_alignment=true
fi
if eval $conditions_to_align; then
  ${scripts_dir}/unifier.py \
    --conversion input \
    --json_params $parameters_file \
    --to_convert $sequence_file \
    --tool AlphaFold3 \
    || { echo "Input conversion failed. Exiting."; exit 1; }
  echo "Running alignment for $sequence_name"
  if $only_msas; then
    following_msas=false
  else
    following_msas=true
  fi
  ${scripts_dir}/create_jobfile.py \
  --job_type=alignment \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file} \
  --mf_following_msas=${following_msas} \
  --tool=${tool} \
  || { echo "Creating alignement jobfile failed. Exiting."; exit 1; }

  ALIGNMENT_ID=$(sbatch --parsable ${sequence_name}_${run_name}_alignment.slurm)
  waiting_for_alignment=true
  if $only_msas; then
    mkdir -p ${logs_dir}/${sequence_name}/${run_name}/
    mv ${sequence_name}_${run_name}_* ${logs_dir}/${sequence_name}/${run_name}/
    echo "Only run sequence alignment."
    exit 1
  fi
elif [[ ! -f $msas_precomputed/msas_alphafold3_data.json ]]; then
    echo "Directory $msas_precomputed does not exits or does not contain msas."
    exit 1
else
  echo "$msas_precomputed are valid."
  mkdir -p ${output_dir}/${sequence_name}/${run_name}
  if ! $waiting_for_alignment; then
    ${scripts_dir}/unifier.py \
      --conversion input_inference \
      --to_convert $msas_precomputed/msas_alphafold3_data.json \
      --json_params $parameters_file \
      --batches_file ${sequence_name}_${run_name}_batches.json \
      --tool AlphaFold3 \
      || { echo "Input preparation for inference has failed. Exiting."; exit 1; }
  fi
fi

# Create and launch inference jobarray
${scripts_dir}/create_jobfile.py \
  --job_type=jobarray \
  --sequence_name=${sequence_name} \
  --run_name=${run_name} \
  --path_to_parameters=${parameters_file} \
  --mf_before_inference $using_jobid \
  --tool AlphaFold3 \
  || { echo "Jobfiles creation has failed. Exiting."; exit 1; }

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
  --path_to_parameters=${parameters_file} \
  --tool AlphaFold3 \
|| { echo "Jobfiles creation has failed. Exiting."; exit 1; }

sbatch --dependency=afterok:$ARRAY_ID ${sequence_name}_${run_name}_post_treatment.slurm

# Store jobiles and batches elements in logs
mkdir -p ${logs_dir}/${sequence_name}/${run_name}/
mv ${sequence_name}_${run_name}_* ${logs_dir}/${sequence_name}/${run_name}/
