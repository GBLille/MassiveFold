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

def lighten_all_pkl(directory):
  try:
    os.mkdir(f'{directory}/light_pkl')
  except OSError as error:
    shutil.rmtree(f'{directory}/light_pkl')
    os.mkdir(f'{directory}/light_pkl')

  to_keep = [
    "num_recycles", "predicted_aligned_error",
    "predicted_lddt", "plddt", "ptm", "iptm",
    "ranking_confidence", "max_predicted_aligned_error"]

  pkl_files = [ file for file in os.listdir(directory) if (file.startswith('result') and file.endswith('.pkl')) ]

  for pkl in pkl_files:
    with open(f"{directory}/{pkl}", 'rb') as pickle_input:
      content = pickle.load(pickle_input)

    lighten_pkl(content, to_keep, pkl, directory)

if __name__ == '__main__':
  try:
    directory = sys.argv[1]
  except:
    print('Usage ./lighten_pkl.py <PATH_TO_OUTPUT>')
    sys.exit()

  if sys.argv[1] == '-h':
    print('Usage ./lighten_pkl.py <PATH_TO_OUTPUT>')
    sys.exit()

  print(f'Extracting pkl from {os.path.abspath(directory)}')
  lighten_all_pkl(directory)
