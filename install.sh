#!/bin/bash

install_path=$1
alphafold_databases=$2

if [ -z $install_path ] || [ ! -d $install_path ]; then
  echo "Path to installation is missing or does not exist"
  echo "$install_path"
  exit 1
fi

if [ -z $alphafold_databases ]; then
  echo "Give a valid path for alphafold databases"
  exit 1
elif [ ! -d $alphafold_databases ]; then
  echo "$alphafold_databases doesn't exists"
  exit 1
fi

mkdir -p $install_path/massivefold_runs/input
mkdir $install_path/massivefold_runs/output
mkdir $install_path/massivefold_runs/log

pipeline=$install_path/massivefold_runs/scripts/
mkdir $pipeline
cp -r MF_scripts/parallelization/* $pipeline

param_file=$pipeline/params.json

params_with_paths=$(cat $param_file | python3 -c "
import json
import sys

params = json.load(sys.stdin)

params['MF_parallel']['run_massivefold'] = '$(realpath ./run_alphafold.py)'
params['MF_parallel']['run_massivefold_plots'] = '$(realpath ./MF_scripts/plots/MF_plots.py)'
params['MF_parallel']['data_dir'] = '$(echo $alphafold_databases)'
params['MF_parallel']['jobfile_headers_dir'] = './headers/'
params['MF_parallel']['jobfile_templates_dir'] = './templates/'
params['MF_parallel']['output_dir'] = '../output/'
params['MF_parallel']['logs_dir'] = '../log/'
params['MF_parallel']['input_dir'] = '../input/'
with open('$param_file', 'w') as params_output:
  json.dump(params, params_output, indent=4)
"
)

cat $param_file
