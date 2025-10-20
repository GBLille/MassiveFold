#!/usr/bin/env python

import math
import json
import numpy as np
import pandas as pd
from copy import deepcopy
from absl import app, flags
import sys

flags.DEFINE_integer('predictions_per_model', 25, 
                     'Choose the number of predictions inferred by each neural network model.')
flags.DEFINE_integer('batch_size', 25, 
                     'Standard size of a prediction batch, if the number of prediction per model\
                       is not a multiple of it, the last batch will be smaller .')
flags.DEFINE_list('models_to_use', None, 'Select the models used for prediction among the five models of each AlphaFold2 version (15 in total).')
flags.DEFINE_string('sequence_name', '', 'Name of the sequence to predict.')
flags.DEFINE_string('run_name', '', 'Name of the run.')
flags.DEFINE_string('path_to_parameters', '', 'Parameters to use, contains models_to_use')
flags.DEFINE_string('to_screen', '', 'Csv file that contains all ligands to screen')
flags.DEFINE_enum('tool', 'AFmassive', ['AFmassive', 'AlphaFold3', 'ColabFold'], 'Specify the tool that you want to use')
FLAGS = flags.FLAGS

def batches_per_model(preds_per_model:int):
  batch_nb = math.ceil(preds_per_model/FLAGS.batch_size)
  batch_sizes = []
  for _ in range(1, batch_nb+1):
    # split total by batch of the same size, the remaining is a single smaller batch
    if preds_per_model - FLAGS.batch_size >= 0:
      preds_per_model -= FLAGS.batch_size
      batch_sizes.append(FLAGS.batch_size)
    else:
      batch_sizes.append(preds_per_model)

  batch_edges = list(np.cumsum([0] + batch_sizes))
  one_model_batches = {
    i+1: {'start': str(batch_edges[i]), 'end': str(batch_edges[i+1]-1)}
    for i in range(len(batch_edges)-1)
  }
  return one_model_batches

def batches_per_ligand(ligands, preds_per_model:int):
  df = pd.read_csv(ligands)
  smiles_and_ccdcode = []

  batches = {}
  for i, (id, smiles, ccdcode) in enumerate(zip(df['id'], df["smiles"], df["ccdCode"])):
    smiles = "" if pd.isna(smiles) else smiles
    ccdcode = "" if pd.isna(ccdcode) else ccdcode
    smiles, ccdcode = smiles.strip(), ccdcode.strip()
    batches[str(i)] = { 'start': 0, 'end': preds_per_model - 1, "id": id.lower(), "smiles": smiles, "ccdcode": [ ccdcode ]}

    if smiles and ccdcode:
      smiles_and_ccdcode.append(i)

  if smiles_and_ccdcode:
    print("For each ligand, put either a smiles or a ccdCodes but not both, you gave both for the following lines:")
    print(", ".join([ str(i) for i in smiles_and_ccdcode ]))
    sys.exit()

  return batches

def batches_all_models(batches_unit, all_models):
  batches = {}
  for i, model in enumerate(all_models):
    unadded_batch = deepcopy(batches_unit)
    for batch in unadded_batch:
      unadded_batch[batch].update({'model': model})
      batches[str(batch+i*len(batches_unit)-1)] = unadded_batch[batch]

  return batches

def main(argv):
  if not FLAGS.path_to_parameters:
    raise ValueError('--path_to_parameters is required')

  with open(FLAGS.path_to_parameters, 'r') as params:
    all_params = json.load(params)

  models = []
  tool_code = ""
  if FLAGS.tool == "AFmassive":
    models = all_params['massivefold']['models_to_use']
    tool_code = "AFM"
  elif FLAGS.tool == "AlphaFold3":
    models = ["AlphaFold3"]
    tool_code = "AF3"
  elif FLAGS.tool == "ColabFold":
    models = all_params['massivefold']['models_to_use']
    tool_code = "CF"
  else:
    print("Tool should be one of the three following: AFmassive, AlphaFold3 or ColabFold.")
    sys.exit()

  if FLAGS.tool == "AlphaFold3":
    model_names = ["AlphaFold3"]

  else:
    model_preset = all_params[f'{tool_code}_run'][f'model_preset']  
    model_names = []
    if model_preset == 'multimer':
      model_names = [
      'model_1_multimer_v1', 'model_2_multimer_v1', 'model_3_multimer_v1', 'model_4_multimer_v1', 'model_5_multimer_v1',
      'model_1_multimer_v2', 'model_2_multimer_v2', 'model_3_multimer_v2', 'model_4_multimer_v2', 'model_5_multimer_v2',
      'model_1_multimer_v3', 'model_2_multimer_v3', 'model_3_multimer_v3', 'model_4_multimer_v3', 'model_5_multimer_v3'
      ]
    elif model_preset == 'monomer_ptm':
      model_names = [ 'model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm', 'model_5_ptm' ]

  if models:
    model_names = [model for model in model_names if model in models]
  if FLAGS.models_to_use:
    model_names = [model for model in model_names if model in FLAGS.models_to_use]

  print(f"Running inference on models: {(', ').join(model_names)}") 
  print(f"Running {FLAGS.predictions_per_model} predictions on each of the {len(model_names)} models")
  print(f"Total prediction number: {FLAGS.predictions_per_model * len(model_names)}")

  if not FLAGS.to_screen:
    # Divide the predictions in batches 
    per_model_batches = batches_per_model(preds_per_model=FLAGS.predictions_per_model)
    # Distribute the batches on all models
    all_model_batches = batches_all_models(per_model_batches, model_names)
  else:
    all_model_batches = batches_per_ligand(FLAGS.to_screen, FLAGS.predictions_per_model)

  with open(f"{FLAGS.sequence_name}_{FLAGS.run_name}_batches.json", "w") as json_output:
    json.dump(all_model_batches, json_output, indent=4)

if __name__ == "__main__":
  app.run(main)
