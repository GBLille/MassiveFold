#!/bin/bash

alphafold_databases=$1

host=$(hostname | cut -c1-8)

if [ "$host" == 'jean-zay' ]; then
  host_is_jeanzay=true
  echo "Currently on Jean Zay cluster, using prebuilt headers and json parameter files."
else
  host_is_jeanzay=false
fi

if ! $host_is_jeanzay; then
  if [ -z $alphafold_databases ]; then
    echo "Give a valid path for alphafold databases"
    exit 1
  elif [ ! -d $alphafold_databases ]; then
    echo "$alphafold_databases doesn't exists"
    exit 1
  fi

  # create massivefold env
  conda="$(conda info --base)/etc/profile.d/conda.sh"
  source $conda
  conda env create -f environment.yml
  conda activate massivefold
  wget -O ${CONDA_PREFIX}/lib/python3.8/site-packages/alphafold/common/stereo_chemical_props.txt https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

  # add run_AFmassive.py and massivefold_plots.py in path (python executables)
  wget -O $CONDA_PREFIX/bin/run_AFmassive.py https://raw.githubusercontent.com/GBLille/AFmassive/main/run_AFmassive.py
  chmod +x $CONDA_PREFIX/bin/run_AFmassive.py

  cp -r massivefold/plots $CONDA_PREFIX/bin/
  cp massivefold/massivefold_plots.py $CONDA_PREFIX/bin/
  chmod +x $CONDA_PREFIX/bin/massivefold_plots.py
fi

# set file tree
runs=massivefold_runs
mkdir -p $runs/input
mkdir $runs/output
mkdir $runs/log
cp examples/H1144.fasta $runs/input

# scripts and files for each pipeline (currently only AFmassive)
cp massivefold/run_massivefold.sh $runs
cp -r massivefold/parallelization/headers $runs

if $host_is_jeanzay; then
  cp massivefold/parallelization/jeanzay_params.json $runs/AFmassive_params.json
  cat $runs/AFmassive_params.json
  echo "Taking Jean Zay's prebuilt header and renaming them."
  mv $runs/headers/example_header_alignment_jeanzay.slurm $runs/headers/alignment.slurm
  mv $runs/headers/example_header_jobarray_jeanzay.slurm $runs/headers/jobarray.slurm
  mv $runs/headers/example_header_post_treatment_jeanzay.slurm $runs/headers/post_treatment.slurm
  exit 1
fi

param_file=$runs/AFmassive_params.json
cp massivefold/parallelization/params.json $param_file

# parameters auto setting
params_with_paths=$(cat $param_file | python3 -c "
import json
import sys

params = json.load(sys.stdin)

params['massivefold']['run_massivefold'] = 'run_AFmassive.py'
params['massivefold']['run_massivefold_plots'] = '$CONDA_PREFIX/bin/massivefold_plots.py'
#params['massivefold']['run_massivefold_plots'] = '$(realpath ./massivefold/massivefold_plots.py)'
params['massivefold']['data_dir'] = '$(realpath $alphafold_databases)'
params['massivefold']['jobfile_templates_dir'] = '../massivefold/parallelization/templates'
params['massivefold']['scripts_dir'] = '../massivefold/parallelization'
params['massivefold']['jobfile_headers_dir'] = './headers'
params['massivefold']['output_dir'] = './output'
params['massivefold']['logs_dir'] = './log'
params['massivefold']['input_dir'] = './input'
with open('$param_file', 'w') as params_output:
  json.dump(params, params_output, indent=4)
"
)

cat $param_file
