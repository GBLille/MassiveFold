#!/usr/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
raw_confidence="https://raw.githubusercontent.com/samuelmurail/alphafold/refs/heads/chain_iptm/alphafold/common/confidence.py"
raw_modules_multimer="https://raw.githubusercontent.com/samuelmurail/alphafold/refs/heads/chain_iptm/alphafold/model/modules_multimer.py"


alphafold_path=${CONDA_PREFIX}/lib/python3.11/site-packages/alphafold

wget -O ${alphafold_path}/common/confidence.py $raw_confidence
wget -O ${alphafold_path}/model/modules_multimer.py $raw_modules_multimer
