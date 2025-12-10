#!/usr/bin/env python

import pickle
import json
import os
import sys
import shutil
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('path_to_output')
parser.add_argument('--pickle_size', 
                    default="full",
                    choices=["full", "light", "custom", "delete"],
                    help='How to treat the stored pickles in the output.',
                    required=True)
parser.add_argument('--parameters', help="Json file containing the parameters for custom pickle size (--pickle_size=custom).")
parser.add_argument('--keep_full_pickles',  action="store_true")

def extract_af3_batch_input_msas(directory: str, json_files: list):

  entity_count = {}
  msas_templates_paths = {}
  total = len(json_files)
  print(f"Now processing {total} af3 batch files in {directory}")
  for file_ind, file in enumerate(json_files):

    filename = os.path.join(directory, file)
    data = json.load(open(filename, 'r'))
    with open(filename, 'r') as batch_file:
      data = json.load(batch_file)

    # extract msas/templates from first json batch and write it to .a3m/mmcif files
    if file_ind == 0:
      entity_count.update({ list(data["sequences"][nseq].keys())[0]: 0 for nseq in range(len(data["sequences"])) })
      for seq_ind, seq in enumerate(data["sequences"]):
        entity = list(seq.keys())[0]
        entity_count[entity] += 1
        msas_templates_paths[f"{entity}_{entity_count[entity]}"] = { }
        msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"] = {}

        # first the msas
        for msas in ["unpairedMsa", "pairedMsa"]:
          # skip chains with no alignments (ligands)
          if not msas in data["sequences"][seq_ind][entity] or not data["sequences"][seq_ind][entity][msas]:
            continue
          fileout = os.path.join(directory, f"{entity}_{entity_count[entity]}_{msas}.a3m")
          msas_templates_paths[f"{entity}_{entity_count[entity]}"][msas] = fileout
          with open(fileout, 'w') as msas_file:
            msas_file.write(data["sequences"][seq_ind][entity][msas])
        # then the templates
        if not "templates" in data["sequences"][seq_ind][entity] or not data["sequences"][seq_ind][entity]["templates"]:
          continue
        templates = data["sequences"][seq_ind][entity]["templates"]
        msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"] = {}
        for n, template in enumerate(templates):
          fileout = os.path.join(directory, f"{entity}_{entity_count[entity]}_template_{n}.mmcif")
          msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"][n] = fileout
          with open(fileout, 'w') as template_file:
            template_file.write(data["sequences"][seq_ind][entity]["templates"][n]["mmcif"])

    # reference the .a3m/mmcif paths in the right entities fields
    entity_count = { list(data["sequences"][nseq].keys())[0]: 0 for nseq in range(len(data["sequences"])) }
    for seq_ind, seq in enumerate(data["sequences"]):
      entity = list(seq.keys())[0]
      entity_count[entity] += 1
      # first the msas
      for msas in ["unpairedMsa", "pairedMsa"]:
        if not msas in data["sequences"][seq_ind][entity] or not data["sequences"][seq_ind][entity][msas]:
          continue
        data["sequences"][seq_ind][entity][msas] = ""
        path_of_msas = msas_templates_paths[f"{entity}_{entity_count[entity]}"][msas]
        data["sequences"][seq_ind][entity][f"{msas}Path"] = path_of_msas

        # then the templates
        if not "templates" in data["sequences"][seq_ind][entity] or not data["sequences"][seq_ind][entity]["templates"]:
          continue
        templates = data["sequences"][seq_ind][entity]["templates"]
        for n, template in enumerate(templates):
          data["sequences"][seq_ind][entity]["templates"][n]["mmcif"] = ""
          path_of_mmcif = msas_templates_paths[f"{entity}_{entity_count[entity]}"]["templates"][n]
          data["sequences"][seq_ind][entity]["templates"][n]["mmcifPath"] = path_of_mmcif

    json.dump(data, open(filename, 'w'), indent=4)

    bar_length = 50
    progress = (file_ind + 1) / total
    filled = int(progress * bar_length)
    bar = "█" * filled + "-" * (bar_length - filled)
    percent = int(progress * 100)
    print(f"\r|{bar}| {percent:3d}% ({file_ind+1}/{total})", end="", flush=True)
  print()

def format_entry(key: str, value, formats):
  possible_formats = [ "npfloat32", "lst" ]
  formattable_entries = [ "predicted_aligned_error", "plddt" ]
  # no format specified for the key
  if key not in formats:
    return value
  if key not in formattable_entries:
    print(f"Key {key} is not formattable, skiping.")
    return value

  user_format = formats[key]
  if user_format not in possible_formats:
    raise ValueError(f"{formats[key]} not a valid format ({', '.join(possible_formats)})")

  formatted_value = value
  if user_format == "npfloat32":
    formatted_value = np.array(value, dtype="float32")
  elif user_format == "lst":
    formatted_value = list(value)

  return formatted_value

def lighten_single_pkl(pkl, directory, parameters):
  with open(f"{directory}/{pkl}", 'rb') as pickle_input:
    initial_content = pickle.load(pickle_input)
  # lighten the pkl content
  content = {}
  for elem in initial_content:
    if elem not in parameters["keys"]:
      continue
    formatted_value = format_entry(elem, initial_content[elem], parameters["format"])
    content[elem] = formatted_value
  # write out the lightened pkl
  with open(f"{directory}/light_pkl/{pkl}", 'wb') as pickle_output:
    pickle.dump(content, pickle_output)
  return content

def delete_pickles(pkl_files, directory):
  total = len(pkl_files)
  print(f"Now deleting {total} pickles in {directory}")
  for i, pkl in enumerate(pkl_files):
    os.remove(os.path.join(directory, pkl))
    # hand-made progress bar
    bar_length = 50
    progress = (i + 1) / total
    filled = int(progress * bar_length)
    bar = "█" * filled + "-" * (bar_length - filled)
    percent = int(progress * 100)
    print(f"\r|{bar}| {percent:3d}% ({i+1}/{total})", end="", flush=True)
  print()

def lighten_all_pkl(directory, parameters, to_json: bool):
  directory_content = os.listdir(directory)

  # delete screening pkls, need a refactor (should be called on each ligand directory instead)
  if os.path.exists(os.path.join(directory, 'screening_inputs')):
    print('AF3 screening detected')
    if args.pickle_size != "delete":
      sys.exit()
    inputs = os.path.join(directory, 'screening_inputs')
    print(f"Remove pickles from directory that inputs in {inputs} have output.")
    ligands_path = os.path.dirname(inputs)

    # remove all ligands pickles
    all_ligands = [ json.load(open(os.path.join(inputs, i), 'r'))['name'] for i in os.listdir(inputs) ]
    for ligand in all_ligands:
      pickles = [ i for i in os.listdir(os.path.join(ligands_path, ligand)) if i.endswith('.pkl') ]
      for pkl in pickles:
        pkl_path = os.path.join(ligands_path, ligand, pkl)
        os.remove(pkl_path)
    sys.exit()

  pkl_files = [ file for file in directory_content if (file.startswith('result') and file.endswith('.pkl')) ]

  if not pkl_files:
    print(f"No pickle files detected in {directory}")

    return

  if os.path.isdir(f'{directory}/light_pkl'):
    shutil.rmtree(f'{directory}/light_pkl')
  os.mkdir(f'{directory}/light_pkl')

  if to_json:
    json_dir = f'{directory}/json_output'
    if os.path.isdir(json_dir):
      shutil.rmtree(json_dir)
    os.mkdir(json_dir)
  """
  try:
    os.mkdir(f'{directory}/light_pkl')
  except OSError as error:
    shutil.rmtree(f'{directory}/light_pkl')
    os.mkdir(f'{directory}/light_pkl')
  """
  total = len(pkl_files)
  for i, pkl in enumerate(pkl_files):
    pkl_content = lighten_single_pkl(pkl, directory, parameters)
    jsonified = {}
    for key in pkl_content:
      if isinstance(pkl_content[key], np.ndarray):
        jsonified[key] = pkl_content[key].tolist()
      else:
        jsonified[key] = pkl_content[key]

    if to_json:
      json.dump(jsonified, open(f"{json_dir}/{pkl.replace('.pkl', '.json')}", 'w'), indent=4)


    # hand-made progress bar
    bar_length = 50
    progress = (i + 1) / total
    filled = int(progress * bar_length)
    bar = "█" * filled + "-" * (bar_length - filled)
    percent = int(progress * 100)
    print(f"\r|{bar}| {percent:3d}% ({i+1}/{total})", end="", flush=True)
  print()

  return pkl_files

if __name__ == '__main__':
  args = parser.parse_args()
  directory = args.path_to_output
  to_json = False # development in progress

  parameters = {
    "keys": [
        "num_recycles", "predicted_aligned_error", "predicted_lddt",
        "plddt", "ptm", "iptm", "ranking_confidence", "max_predicted_aligned_error"
      ],
    "format": {
      "predicted_aligned_error": "npfloat32", # npfloat32 | lst
      "plddt": "npfloat32"
    }
  }

  af3_batch_files = [ i for i in os.listdir(directory) if i.startswith('af3') and i.endswith('.json') ]
  if af3_batch_files:
    print("Detected AlphaFold3 output.")
    print("Reference large input elements out of the af3 batches files")
    extract_af3_batch_input_msas(directory, af3_batch_files)

  # step to lighten (or not) the pickles
  if args.pickle_size == "full":
    print("No modification of the pickle files")
    if to_json:
      json_dir = f'{directory}/json_output'
      if os.path.isdir(json_dir):
        shutil.rmtree(json_dir)
      os.mkdir(json_dir)
      print(f"Converting the full pickle files to json at {json_dir}")
      pkl_files = [ file for file in os.listdir(directory) if (file.startswith('result') and file.endswith('.pkl')) ]
      total = len(pkl_files)
      for i, pkl in enumerate(pkl_files):
        pkl_content = pickle.load(open(f"{directory}/{pkl}", 'rb'))
        jsonified = {}
        for key in pkl_content:
          if isinstance(pkl_content[key], np.ndarray):
            jsonified[key] = pkl_content[key].tolist()
          elif isinstance(pkl_content[key], np.float64):
            jsonified[key] = float(pkl_content[key])
          else:
            jsonified[key] = pkl_content[key]
            print(jsonified[key])

        json.dump(jsonified, open(f"{json_dir}/{pkl.replace('.pkl', '.json')}", 'w'), indent=4)
        bar_length = 50
        progress = (i + 1) / total
        filled = int(progress * bar_length)
        bar = "█" * filled + "-" * (bar_length - filled)
        percent = int(progress * 100)
        print(f"\r|{bar}| {percent:3d}% ({i+1}/{total})", end="", flush=True)
      print()

  elif args.pickle_size in [ "light", "custom" ]:
    delete_after_lightening = not args.keep_full_pickles
    print(f'Extracting pkl from {os.path.abspath(directory)}')
    # read user's custom parameters for the lighter pickles
    if args.pickle_size == "custom":
      parameters = json.load(open(args.parameters, "r"))
      print(f"Using custom parameters:\n{parameters}")
    pickles = lighten_all_pkl(directory, parameters, to_json=False)
    if delete_after_lightening and pickles != None:
      delete_pickles(pickles, directory)

  elif args.pickle_size == "delete":
    pickles = [ pkl for pkl in os.listdir(directory) if pkl.endswith('.pkl') ]
    delete_pickles(pickles, directory)
