#!/usr/bin/env python

import argparse
import math
import json
import numpy as np
import pandas as pd
from copy import deepcopy
import sys
import os

def parse_csv(value):
  if value is None:
    return None
  if value == '':
    return []
  return [item.strip() for item in value.split(',') if item.strip()]

def batches_per_model(preds_per_model, batch_size):
  batch_nb = math.ceil(preds_per_model/batch_size)
  batch_sizes = []
  for _ in range(1, batch_nb+1):
    # split total by batch of the same size, the remaining is a single smaller batch
    if preds_per_model - batch_size >= 0:
      preds_per_model -= batch_size
      batch_sizes.append(batch_size)
    else:
      batch_sizes.append(preds_per_model)

  batch_edges = list(np.cumsum([0] + batch_sizes))
  one_model_batches = {
    i+1: {'start': str(batch_edges[i]), 'end': str(batch_edges[i+1]-1)}
    for i in range(len(batch_edges)-1)
  }
  return one_model_batches

def batches_per_ligand(ligands, preds_per_model):
  df = pd.read_csv(ligands)
  multiple_types = []

  batches = {}
  for i, (id, smiles, ccdcode, iupac) in enumerate(zip(df['id'], df["smiles"], df["ccdCode"], df["IUPAC"])):
    smiles = "" if pd.isna(smiles) else smiles
    ccdcode = "" if pd.isna(ccdcode) else ccdcode
    iupac = "" if pd.isna(iupac) else iupac
    smiles, ccdcode, iupac = smiles.strip(), ccdcode.strip(), iupac.strip()
    batches[str(i)] = {
      'start': 0, 'end': preds_per_model - 1, "id": id.lower(),
      "smiles": smiles, "ccdcode": [ ccdcode ], "iupac": iupac
    }

    if sum([bool(smiles), bool(ccdcode), bool(iupac)]) > 1:
      multiple_types.append(i)

  if multiple_types:
    print("For each ligand, put either a smiles, ccdCodes or IUPAC but not more than one, you gave multiple for the following lines:")
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

def detect_model_preset(fasta_file):
  with open(fasta_file, 'r') as fasta_content:
    lines = fasta_content.readlines()

  model_preset = "monomer_ptm"
  sequence_number = 0
  for line in lines:
    if line.startswith('>'):
      sequence_number += 1
    if sequence_number > 1:
      model_preset = "multimer"
      break
  return model_preset

def main(
  predictions_per_model,
  batch_size,
  models_to_use,
  sequence_name,
  run_name,
  path_to_parameters,
  to_screen,
  tool):
  if not path_to_parameters:
    raise ValueError('--path_to_parameters is required')

  with open(path_to_parameters, 'r') as params:
    all_params = json.load(params)

  model_preset =  detect_model_preset(
    os.path.join(all_params["massivefold"]["input_dir"], f"{sequence_name}.fasta")
  )

  tool_code = ""
  if tool == "AFmassive":
    tool_code = "AFM"
  elif tool == "AlphaFold3":
    models = ["AlphaFold3"]
    tool_code = "AF3"
  elif tool == "ColabFold":
    tool_code = "CF"
  else:
    print("Tool should be one of the three following: AFmassive, AlphaFold3 or ColabFold.")
    sys.exit()

  if tool == "AlphaFold3":
    model_names = ["AlphaFold3"]

  else:
    models = []
    models_string = all_params['massivefold']["models_to_use"]
    models = models_string.split(',') if models_string else []
    model_names = []
    if model_preset == 'multimer':
      model_names = [
      'model_1_multimer_v1', 'model_2_multimer_v1', 'model_3_multimer_v1', 'model_4_multimer_v1', 'model_5_multimer_v1',
      'model_1_multimer_v2', 'model_2_multimer_v2', 'model_3_multimer_v2', 'model_4_multimer_v2', 'model_5_multimer_v2',
      'model_1_multimer_v3', 'model_2_multimer_v3', 'model_3_multimer_v3', 'model_4_multimer_v3', 'model_5_multimer_v3'
      ]
    elif model_preset == 'monomer_ptm':
      model_names = [ 'model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm', 'model_5_ptm' ]

    non_existing_models = []
    non_existing_models.extend([ i for i in models if i not in model_names ])
    non_existing_models.extend([ i for i in (models_to_use or []) if i not in model_names ])

    if non_existing_models:
      raise ValueError(f"Model '{', '.join(non_existing_models)}' does not exist for preset '{model_preset}'")

    if models:
      model_names = [model for model in model_names if model in models]
    if models_to_use:
      model_names = [model for model in model_names if model in models_to_use]

  print(f"Running inference on models: {(', ').join(model_names)}") 
  print(f"Running {predictions_per_model} predictions on each of the {len(model_names)} models")
  print(f"Total prediction number: {predictions_per_model * len(model_names)}")

  if not to_screen:
    # Divide the predictions in batches 
    per_model_batches = batches_per_model(preds_per_model=predictions_per_model, batch_size=batch_size)
    # Distribute the batches on all models
    all_model_batches = batches_all_models(per_model_batches, model_names)
  else:
    all_model_batches = batches_per_ligand(to_screen, predictions_per_model)

  with open(f"{sequence_name}_{run_name}_batches.json", "w") as json_output:
    json.dump(all_model_batches, json_output, indent=4)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--predictions_per_model',
    type=int,
    default=25,
    help='Choose the number of predictions inferred by each neural network model.')
  parser.add_argument(
    '--batch_size',
    type=int,
    default=25,
    help='Standard size of a prediction batch, if the number of prediction per model is not a multiple of it, the last batch will be smaller.')
  parser.add_argument(
    '--models_to_use',
    type=parse_csv,
    default=None,
    help='Select the models used for prediction among the five models of each AlphaFold2 version (15 in total).')
  parser.add_argument('--sequence_name', default='', help='Name of the sequence to predict.')
  parser.add_argument('--run_name', default='', help='Name of the run.')
  parser.add_argument('--path_to_parameters', default='', help='Parameters to use, contains models_to_use')
  parser.add_argument('--to_screen', default='', help='Csv file that contains all ligands to screen')
  parser.add_argument(
    '--tool',
    default='AFmassive',
    choices=['AFmassive', 'AlphaFold3', 'ColabFold'],
    help='Specify the tool that you want to use')

  parsed = parser.parse_args()
  main(
    parsed.predictions_per_model,
    parsed.batch_size,
    parsed.models_to_use,
    parsed.sequence_name,
    parsed.run_name,
    parsed.path_to_parameters,
    parsed.to_screen,
    parsed.tool
  )
