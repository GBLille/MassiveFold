#!/bin/bash

install_env () {
  env=$1
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if [[ $env == "massivefold" ]]; then
    echo "Installing MassiveFold environment"
    conda env create -f environment.yml
    conda activate massivefold
    conda config --env --set channel_priority flexible
    python -m pip install -e .
  elif [[ $env == "colabfold" ]]; then
    echo "Installing ColabFold environment"
    conda activate massivefold || { echo "massivefold environment is needed and is missing"; exit 1; }
    CONDA_OVERRIDE_CUDA="11.8" conda env create -f mf-colabfold.yml
  elif [[ $env == "afmassive" ]]; then
    echo "Installing afmassive environment"
    conda activate massivefold || { echo "massivefold environment is needed and is missing"; exit 1; }
    CONDA_OVERRIDE_CUDA="11.8" conda env create -f mf-afmassive.yml
    conda activate mf-afmassive-1.1.6
    wget -O "${CONDA_PREFIX}/lib/python3.10/site-packages/alphafold/common/stereo_chemical_props.txt" https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
    wget -O "${CONDA_PREFIX}/bin/run_AFmassive.py" https://raw.githubusercontent.com/GBLille/AFmassive/v1.1.6/run_AFmassive.py
    chmod +x "${CONDA_PREFIX}/bin/run_AFmassive.py"
    conda deactivate
  elif [[ $env == "alphafold3" ]]; then
    echo "Installing alphafold3 environment"
    conda activate massivefold || { echo "massivefold environment is needed and is missing"; exit 1; }
    conda env create -f mf-alphafold3.yml
    conda activate mf-alphafold-3.0.1
    build_data
    wget -O "${CONDA_PREFIX}/bin/run_alphafold.py" https://raw.githubusercontent.com/google-deepmind/alphafold3/v3.0.1/run_alphafold.py
    sed -i '1i #!/usr/bin/env python' "${CONDA_PREFIX}/bin/run_alphafold.py"
    chmod +x "${CONDA_PREFIX}/bin/run_alphafold.py"
    conda deactivate
  fi
}

do_help=false
db_af=false
db_cf=false
db_af3=false
only_create_env=false
do_not_create_env=false
install_path="massivefold_runs"

while true; do
  case "$1" in
    --alphafold-db)
      alphafold_databases=$2
      db_af=true
      shift 2
      ;;
    --alphafold3-db)
      alphafold3_databases=$2
      db_af3=true
      shift 2
      ;;
    --colabfold-db)
      colabfold_databases=$2
      db_cf=true
      shift 2
      ;;
    --install-path)
      install_path=$2
      shift 2
      ;;
    --no-env)
      do_not_create_env=true
      shift 1
      ;;
    --only-envs)
      only_create_env=true
      shift 1
      ;;
    -h|--help)
      do_help=true
      shift 1
      ;;
    *)
      break
      ;;
  esac
done

USAGE="\
On Jean Zay cluster:\n\
  ./install.sh [--install-path str]\n\
Otherwise:\n\
  ./install.sh [--only-envs] || --alphafold-db str --alphafold3-db str --colabfold-db str [--no-env] [--install-path str]\n\n\
./install.sh -h for more details"

if [[ $do_help == "true" ]]; then
  echo -e "\
Usage:
------
$USAGE
  Options:
    --alphafold-db <str>: path to AlphaFold2 database
    --alphafold3-db <str>: path to AlphaFold3 database
    --colabfold-db <str>: path to ColabFold database
    --install-path <str>: where to create the MassiveFold file architecture (default: massivefold_runs)
    --no-env: do not install the environments, only files and parameters.
      At least one of --alphafold-db or --colabfold-db is required with this option.
    --only-envs: only install the environments (other arguments are not used)"
  exit 1
fi

host=$(hostname | cut -c1-8)
host_is_jeanzay=false
if [[ "$host" == "jean-zay" ]]; then
  host_is_jeanzay=true
fi

if [[ $only_create_env == "true" ]]; then
  install_env "massivefold"
  install_env "alphafold3"
  install_env "afmassive"
  install_env "colabfold"
  echo "Both environments installed, exiting."
  exit 1
fi

if [[ $host_is_jeanzay == "false" ]]; then
  if [[ $db_af == "true" && ! -d $alphafold_databases ]]; then
    echo "$alphafold_databases doesn't exists"
    exit 1
  elif [[ $db_af3 == "true" && ! -d $alphafold3_databases ]]; then
    echo "$alphafold3_databases doesn't exists"
    exit 1
  elif [[ $db_cf == "true" && ! -d $colabfold_databases ]]; then
    echo "$colabfold_databases doesn't exists"
    exit 1
  elif [[ $db_cf == "false" && $db_af == "false" && $db_af3 == "false" ]]; then
    echo -e "$USAGE"
    exit 1
  fi
else
  echo "Currently on Jean Zay cluster, using prebuilt headers and json parameter files."
  mkdir -p "$HOME/af3_datadir/"
  if [[ ! -e "$HOME/af3_datadir/Alphafold3" ]]; then
    ln -s "$DSDIR/Alphafold3/" "$HOME/af3_datadir/"
  fi
fi

conda_sh="$(conda info --base)/etc/profile.d/conda.sh"
source "$conda_sh"

if [[ $do_not_create_env == "false" ]]; then
  install_env "massivefold"
  if [[ $db_af == "true" ]]; then
    install_env "afmassive"
  fi
  if [[ $db_af3 == "true" ]]; then
    install_env "alphafold3"
  fi
  if [[ $db_cf == "true" ]]; then
    install_env "colabfold"
  fi
else
  echo "No env asked, install skipped"
fi


install_cmd=(massivefold install --install-path "$install_path")
if [[ $db_af == "true" ]]; then
  install_cmd+=(--alphafold-db "$alphafold_databases")
fi
if [[ $db_af3 == "true" ]]; then
  install_cmd+=(--alphafold3-db "$alphafold3_databases")
fi
if [[ $db_cf == "true" ]]; then
  install_cmd+=(--colabfold-db "$colabfold_databases")
fi
if [[ $do_not_create_env == "true" ]]; then
  install_cmd+=(--no-env)
fi

"${install_cmd[@]}"
