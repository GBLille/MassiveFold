#!/bin/bash

setup_params () {
  tool=$1
  echo "$tool"
}

install_env () {
  env=$1
  source $(conda info --base)/etc/profile.d/conda.sh
  if [[ $env == "massivefold" ]]; then
    echo "Installing MassiveFold environment"
    CONDA_OVERRIDE_CUDA="11.8" conda env create -f environment.yml
    
    conda activate massivefold
    wget -O ${CONDA_PREFIX}/lib/python3.10/site-packages/alphafold/common/stereo_chemical_props.txt https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

    # add run_AFmassive.py and massivefold_plots.py in path (python executables)
    wget -O $CONDA_PREFIX/bin/run_AFmassive.py https://raw.githubusercontent.com/GBLille/AFmassive/main/run_AFmassive.py
    chmod +x $CONDA_PREFIX/bin/run_AFmassive.py

    cp -r massivefold/plots $CONDA_PREFIX/bin/
    cp massivefold/massivefold_plots.py $CONDA_PREFIX/bin/
    chmod +x $CONDA_PREFIX/bin/massivefold_plots.py

  elif [[ $env == "colabfold" ]]; then
    echo "Installing ColabFold environment"
    CONDA_OVERRIDE_CUDA="11.8" conda env create -f mf_colabfold.yml
  fi
}

# argument parser
while true; do
  case "$1" in
    --alphafold-db)
      alphafold_databases=$2
      db_af=true
      shift 2
      ;;
    --colabfold-db)
      colabfold_databases=$2
      db_cf=true
      shift 2
      ;;
    --no-env)
      do_not_create_env=true
      shift 1
      ;;
    --only-env)
      only_create_env=true
      shift 1
      ;;
    *)
      break
      ;;
  esac
done

host=$(hostname | cut -c1-8)

if $only_create_env; then
  install_env "massivefold"
  install_env "colabfold"
  
  echo "Both environments installed, exiting."
  exit 1
fi

if [ "$host" == 'jean-zay' ]; then
  host_is_jeanzay=true
  echo "Currently on Jean Zay cluster, using prebuilt headers and json parameter files."
else
  host_is_jeanzay=false
fi


if ! $host_is_jeanzay; then
  if [[ $db_af && ! -d $alphafold_databases ]]; then
    echo "$alphafold_databases doesn't exists"
    exit 1
  elif [[ $db_cf && ! -d $colabfold_databases ]]; then
    echo "$colabfold_databases doesn't exists"
    exit 1
  elif [[ ! $db_cf && ! $db_af ]]; then
    echo "Use one of both or both --alphafold-db and --colabfold-db"
    exit 1
  fi
  
  # Install envs only if they are not yet installed
  conda="$(conda info --base)/etc/profile.d/conda.sh"
  source $conda
  mf_env=$(conda env list | grep massivefold | wc -l)
  mf_cf_env=$(conda env list | grep mf-colabfold | wc -l)

  if (( mf_env > 0 )); then  
    echo "massivefold env already installed, skipping this step."
  elif [[ ! $do_not_create_env ]]; then  
    install_env "massivefold"
  fi

  if [[ $db_cf ]] && (( mf_cf_env > 0 )); then
    echo "mf-colabfold env already installed, skipped"
  elif [[ $db_cf ]] || [[ ! $do_not_create_env ]]; then
    install_env "colabfold"
  elif $do_not_create_env; then
    echo "No env asked, install skipped"
  else
    echo "ColabFold databases not provided, install skipped"
  fi
fi

: '
if $db_af; then
  setup_params "afmassive"
fi

if $db_cv; then
  setup_params "colabfold"
fi
'

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

