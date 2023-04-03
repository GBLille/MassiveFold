#!/usr/bin/env python
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Relax AlphaFold protein structure prediction script."""
import json
import os
import re
#import pathlib
import pickle
#import random
#import shutil
import sys
import time
import bz2
#from typing import Dict, Union, Optional

#from absl import app
#from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
#from alphafold.data import pipeline
#from alphafold.data import pipeline_multimer
#from alphafold.data import templates
#from alphafold.data.tools import hhsearch
#from alphafold.data.tools import hmmsearch
from alphafold.model import config
from alphafold.model import model
from alphafold.relax import relax
import numpy as np

from alphafold.model import data
# Internal import (7716).

logging.set_verbosity(logging.INFO)

os.environ['OPENMM_CPU_THREADS']='16'




MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3


def main():
    amber_relaxer = relax.AmberRelaxation(
    max_iterations=RELAX_MAX_ITERATIONS,
    tolerance=RELAX_ENERGY_TOLERANCE,
    stiffness=RELAX_STIFFNESS,
    exclude_residues=RELAX_EXCLUDE_RESIDUES,
    max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS)


    ranking_confidences={}
    unrelaxed_pdbs={}
    relaxed_pdbs={}
    timings={}
    # Load the model outputs.
    result_pickle = sys.argv[1]
    output_dir=os.path.dirname(result_pickle)
    model_name=os.path.basename(result_pickle).replace('result_','').replace('.pkl','')
    relaxed_output_path = result_pickle.replace('result_','relaxed_').replace('.pkl','.pdb')
    if os.path.exists(relaxed_output_path):
        print(f'{relaxed_output_path} exists, delete if you want to rerun.')
        sys.exit(0)
    match=re.search('result_(model_[\d]*["_ptm]*["_multimer"]*["_v2"]*)_\d+\.pkl',result_pickle)
    if match:
        model_name=match.group(1)
    else:
        print(f'Cannot get model_name from {result_pickle}')
        sys.exit()
    
    print(model_name)
    #sys.exit()
    print(os.path.dirname(os.path.realpath(__file__)))
    data_dir=f'{os.path.dirname(os.path.realpath(__file__))}/alphafold_data/'
    print(data_dir)
   
    #model should not be needed but the process_features is a method in the RunModel class
    #the processed_features are used to make a protein class.
    #All of these dependencies could possible be removed, but I wanted the relax to be identical to what is in AF
    #and this was the easiest way
    model_config = config.model_config(model_name)
    model_params = data.get_model_haiku_params(
        model_name=model_name, data_dir=data_dir)
    model_runner = model.RunModel(model_config, model_params)
    
    feature_pickle=f'{output_dir}/features.pkl'
    try:
        f=open(feature_pickle,'rb')
    except:
        print(f'Could not open {feature_pickle}')


    try:
        f=bz2.open(feature_pickle+'.bz2','rb')
    except:
        print(f'Could not open {feature_pickle}.bz2')
        sys.exit()

    processed_feature_dict={}
    feature_dict=pickle.load(f)
    processed_feature_dict = model_runner.process_features(feature_dict, random_seed=42) #mo
    f.close()

    #os.path.join(output_dir, f'result_{model_name}.pkl')

#    try:
#        with open(feature_pickle,'rb') as f:
            
#            processed_feature_dict = model_runner.process_features(feature_dict, random_seed=42) #model_random_seed)
#    except:
        
        
        

    #print(processed_feature_dict)
    with open(result_pickle, 'rb') as f:
      prediction_result=pickle.load(f)
      plddt = prediction_result['plddt']
      ranking_confidences[model_name] = prediction_result['ranking_confidence']
      # Add the predicted LDDT in the b-factor column.
      # Note that higher predicted LDDT value means higher model confidence.
      plddt_b_factors = np.repeat(
        plddt[:, None], residue_constants.atom_type_num, axis=-1)
      unrelaxed_protein = protein.from_prediction(
        features=processed_feature_dict,
        result=prediction_result,
        b_factors=plddt_b_factors,
        remove_leading_feature_dimension=not model_runner.multimer_mode)

      # unrelaxed_pdbs[model_name] = protein.to_pdb(unrelaxed_protein)
      #      unrelaxed_pdb_path = os.path.join(output_dir, f'unrelaxed_{model_name}.pdb')
      #      with open(unrelaxed_pdb_path, 'w') as f:
      #        f.write(unrelaxed_pdbs[model_name])
      #print(unrelaxed_pdbs[model_name])
      #      sys.exit()


#      relaxed_output_path = os.path.join(
#          output_dir, f'relaxed_{model_name}.pdb')
      print(relaxed_output_path)
      #sys.exit()
      t_0 = time.time()
      relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
      timings[f'relax_{model_name}'] = time.time() - t_0
      logging.info(f'Relax took %f s', timings[f'relax_{model_name}'])
      relaxed_pdbs[model_name] = relaxed_pdb_str
      # Save the relaxed PDB.
      with open(relaxed_output_path, 'w') as f:
        f.write(relaxed_pdb_str)



if __name__ == '__main__':
  main()
