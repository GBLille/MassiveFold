#!/usr/bin/env python

import pickle
import json
import os
import shutil
import sys

def lighten_pkl(content, to_keep, file, directory):
  print(file.replace('result_', '').replace('.pkl', ''))
  light_dict = {}

  for i in content:
    if i in to_keep:
      light_dict[i] = content[i]

  with open(f"{directory}/light_pkl/{file}", 'wb') as pickle_output:
    pickle.dump(light_dict, pickle_output)

  return light_dict

def extract_af3_batch_input_msas(directory: str, json_files: list):

  entity_count = {}
  msas_templates_paths = {}
  for i, file in enumerate(json_files):
    print(f"Processing {os.path.join(directory, file)}")
    filename = os.path.join(directory, file)
    data = json.load(open(filename, 'r'))
    with open(filename, 'r') as batch_file:
      data = json.load(batch_file)

    # extract msas/templates from first json batch and write it to .a3m/mmcif files
    if i == 0:
      entity_count.update({ list(data["sequences"][i].keys())[0]: 0 for i in range(len(data["sequences"])) })
      for i, seq in enumerate(data["sequences"]):
        entity = list(seq.keys())[0]
        entity_count[entity] += 1
        msas_templates_paths[f"{entity}_{entity_count[entity]}"] = { }
        msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"] = {}

        # first the msas
        for msas in ["unpairedMsa", "pairedMsa"]:
          # skip chains with no alignments (ligands)
          if not msas in data["sequences"][i][entity] or not data["sequences"][i][entity][msas]:
            continue
          """
          alignments = data["sequences"][i][entity][msas]
          if not alignments:
            continue
          """
          fileout = os.path.join(directory, f"{entity}_{entity_count[entity]}_{msas}.a3m")
          msas_templates_paths[f"{entity}_{entity_count[entity]}"][msas] = fileout
          with open(fileout, 'w') as msas_file:
            msas_file.write(data["sequences"][i][entity][msas])
        # then the templates
        if not "templates" in data["sequences"][i][entity] or not data["sequences"][i][entity]["templates"]:
          continue
        templates = data["sequences"][i][entity]["templates"]
        msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"] = {}
        for n, template in enumerate(templates):
          fileout = os.path.join(directory, f"{entity}_{entity_count[entity]}_template_{n}.mmcif")
          msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"][n] = fileout
          with open(fileout, 'w') as template_file:
            template_file.write(data["sequences"][i][entity]["templates"][n]["mmcif"])

    # reference the .a3m/mmcif paths in the right entities fields
    entity_count = { list(data["sequences"][i].keys())[0]: 0 for i in range(len(data["sequences"])) }
    for i, seq in enumerate(data["sequences"]):
      entity = list(seq.keys())[0]
      entity_count[entity] += 1
      # first the msas
      for msas in ["unpairedMsa", "pairedMsa"]:
        if not msas in data["sequences"][i][entity] or not data["sequences"][i][entity][msas]:
          continue
        data["sequences"][i][entity][msas] = ""
        path_of_msas = msas_templates_paths[f"{entity}_{entity_count[entity]}"][msas]
        data["sequences"][i][entity][f"{msas}Path"] = path_of_msas

        # then the templates
        if not "templates" in data["sequences"][i][entity] or not data["sequences"][i][entity]["templates"]:
          continue
        templates = data["sequences"][i][entity]["templates"]
        for n, template in enumerate(templates):
          data["sequences"][i][entity]["templates"][n]["mmcif"] = ""
          path_of_mmcif = msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"][n]
          data["sequences"][i][entity]["templates"][n]["mmcifPath"] = path_of_mmcif

    json.dump(data, open(filename, 'w'), indent=4)

def lighten_all_pkl(directory):
  directory_content = os.listdir(directory)
  af3_batch_files = [ i for i in directory_content if i.startswith('af3') and i.endswith('.json') ]
  if af3_batch_files:
    extract_af3_batch_input_msas(directory, af3_batch_files)
    sys.exit()
  try:
    os.mkdir(f'{directory}/light_pkl')
  except OSError as error:
    shutil.rmtree(f'{directory}/light_pkl')
    os.mkdir(f'{directory}/light_pkl')

  to_keep = [
    "num_recycles", "predicted_aligned_error",
    "predicted_lddt", "plddt", "ptm", "iptm",
    "ranking_confidence", "max_predicted_aligned_error"]

  pkl_files = [ file for file in directory_content if (file.startswith('result') and file.endswith('.pkl')) ]

  for pkl in pkl_files:
    with open(f"{directory}/{pkl}", 'rb') as pickle_input:
      content = pickle.load(pickle_input)

    lighten_pkl(content, to_keep, pkl, directory)

if __name__ == '__main__':
  if len(sys.argv) != 2:
    print('Usage ./lighten_pkl.py <PATH_TO_OUTPUT>')
    sys.exit()
  directory = sys.argv[1]

  if sys.argv[1] == '-h':
    print('Usage ./lighten_pkl.py <PATH_TO_OUTPUT>')
    sys.exit()

  print(f'Extracting pkl from {os.path.abspath(directory)}')
  lighten_all_pkl(directory)
