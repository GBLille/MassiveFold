#!/usr/bin/env python

import os
import json
import sys
import argparse
import numpy as np
import pandas as pd
from os import symlink
from shutil import copy  as cp, rmtree as rm

parser = argparse.ArgumentParser()
parser.add_argument('--runs_path', help='Path to the runs you want to gather', required=True)
parser.add_argument('--ignore', dest='runs_to_ignore', nargs='+', help="List of run names or" 
"subdirectory in the run path set in --run_path. All these runs (separated by a single space)" 
" won't be taken into account in the runs gathering", default=[])
parser.add_argument('--output_path', help="Path to the output of the runs gathering (default: <RUNS_PATH>/all_pdbs)", required=False)
parser.add_argument('--only_ranking', help="Skips the run gathering, only output a csv ranking file", required=False, action='store_true')
parser.add_argument('--include_pickle', help="If specified, include the .pkl pickle" 
" files in the the gathered results stored in <OUTPUT_PATH>", required=False, action='store_true')

def create_global_ranking(all_runs_path, runs, output_path, ranking_type="debug"):
  map_pred_run = {}
  whole_ranking = {}
  runs = [ os.path.join(all_runs_path, run) for run in runs ]
  ranking_key_score = ''
  try:
    os.mkdir(output_path)
  except FileExistsError:
    pass
  for run in runs:
    ranking_path = os.path.join(run, f'ranking_{ranking_type}.json')
    run_name = os.path.basename(run)
    with open(ranking_path, 'r') as local_ranking_file:
      local_rank = json.load(local_ranking_file)
      if not ranking_key_score:
        ranking_key_score = list(local_rank.keys())[0]
    local_rank[ranking_key_score] = \
    {f"{run_name}_{pred}": local_rank[ranking_key_score][pred] for pred in local_rank[ranking_key_score]}
    local_rank['order'] = [f"{run_name}_{pred}" for pred in local_rank['order']]
    map_pred_run.update({pred: run_name for pred in local_rank['order']})
    whole_ranking.update(local_rank[ranking_key_score])
  sorted_predictions = sorted(whole_ranking.items(), key=lambda x:x[1], reverse=True)
  iptm_ptm = dict(sorted_predictions)
  order = sorted(whole_ranking, reverse=True, key=whole_ranking.get)
  global_ranking = {ranking_key_score: iptm_ptm, 'order': order}
  with open(f"{output_path}/ranking_{ranking_type}.json", 'w') as fileout:
    fileout.write(json.dumps(global_ranking, indent=4)) 
  return map_pred_run

def rank_all(all_runs_path, runs, output_path, ranking_type="debug"):
  runs = [ os.path.join(all_runs_path, run) for run in runs ]
  ranking_key_score = ''
  all_models = pd.DataFrame()
  try:
    os.mkdir(output_path)
  except FileExistsError:
    pass

  for run in runs:
    single_run_models = pd.DataFrame()
    ranking_path = os.path.join(run, f'ranking_{ranking_type}.json')
    with open(ranking_path, 'r') as local_ranking_file:
      local_rank = json.load(local_ranking_file)
    if not ranking_key_score:
      ranking_key_score = list(local_rank.keys())[0]
    
    metric = ranking_key_score
    scores = list(local_rank[ranking_key_score].values())
    model_names = local_rank[ranking_key_score].keys()
    reconstruct_prediction = lambda x, y: f"ranked_{y}_unrelaxed_{x}.pdb"
    predictions = list(map(reconstruct_prediction, model_names, range(len(model_names))))
    
    check_files_existence = lambda x: os.path.isfile(os.path.join(run, x))
  
    files_existence = list(map(check_files_existence, predictions))
    if not all(files_existence):
      print(f'/!\\ some predictions are not found in the run {run}, needs investigation')
      indices = np.logical_not(files_existence).astype(int)
      print(f"Predictions not present: {', '.join(list(np.array(predictions)[indices.astype(bool)]))}\n")
    single_run_models['file'] = predictions
    single_run_models[ranking_key_score] = scores
    parameter_set = os.path.basename(os.path.normpath(run))
    single_run_models['parameters'] = parameter_set
    columns_order = ['parameters', 'file', ranking_key_score]
    single_run_models = single_run_models[columns_order]
    all_models = pd.concat([all_models, single_run_models], axis=0) 
  
  all_models = all_models.sort_values(ranking_key_score, ascending=False, ignore_index=True)
  all_models.to_csv(os.path.join(all_runs_path, 'ranking.csv'), index=None)

  print(all_models)

def move_and_rename(all_runs_path, run_names, output_path, do_include_pickle):
  with open(os.path.join(output_path, 'ranking_debug.json'), 'r') as rank_file:
    global_rank_order = json.load(rank_file)['order']
  for i, prediction in enumerate(global_rank_order):
    # copy the features
    if i == 0:
      features_old_name = os.path.join(all_runs_path, run_names[prediction], "features.pkl")
      features_new_name = os.path.join(output_path, "features.pkl")
      try:
        cp(features_old_name, features_new_name)
      except FileNotFoundError:
        print(f'features.pkl not found in {run_names[prediction]} run')

    # copy the predictions
    pred_first_name = f"unrelaxed_model_{prediction.split('model_')[1]}"
    if os.path.exists(os.path.join(all_runs_path, run_names[prediction], pred_first_name)):
      old_name = f"relaxed_{pred_first_name}.pdb"
    else:
      old_name = f"{pred_first_name}.pdb"

    # Move pdb files and rename with rank
    old_pdb_path = os.path.join(all_runs_path, run_names[prediction], old_name) 
    new_pdb_path = os.path.join(output_path, f"ranked_{i}_{prediction}.pdb")
    cp(old_pdb_path, new_pdb_path)
    
    if do_include_pickle:
      pickles_path = os.path.join(all_runs_path, run_names[prediction])
      local_pred_name = f"model_{prediction.split('model_')[1]}" 
      if os.path.exists(os.path.join(pickles_path, 'light_pkl')):
        pickles_path = os.path.join(pickles_path, 'light_pkl')

      old_pdb_path = os.path.join(pickles_path, f"result_{local_pred_name}.pkl") 
      new_pdb_path = os.path.join(output_path, f"ranked_{i}_{prediction}.pkl")
      cp(old_pdb_path, new_pdb_path)

def create_symlink_without_ranked(all_runs_path, runs):
  for run in runs:
    path = os.path.join(all_runs_path, run)
    ranked_preds = [pred for pred in os.listdir(path) if pred.startswith('ranked_')]
    for pred in ranked_preds:
      old_name = os.path.realpath(os.path.join(path, pred))
      new_name = os.path.join(path, pred.split('_', 2)[-1])
      try:
        os.symlink(old_name, new_name)
      except FileExistsError:
        pass

def delete_symlinks(all_runs_path, runs):
  for run in runs:
    path_run = os.path.join(all_runs_path, run)
    files = os.listdir(path_run)
    for symlink in files:
      symlink_file = os.path.join(path_run, symlink)
      if os.path.islink(symlink_file):
        os.remove(symlink_file)

def check_all_runs(all_runs_path, ignored_directories, ranking_type='debug'):
  considered_runs = []
  
  for run in os.listdir(all_runs_path):
    if run not in ignored_directories:
      ranking_path = os.path.join(all_runs_path, run, f'ranking_{ranking_type}.json')
      if os.path.isfile(ranking_path):
        try:
          with open(ranking_path, 'r') as ranking_file:
            json.load(ranking_file)
          considered_runs.append(run)
        except:
          print(f'Something went wrong with ranking_debug.json for run {run}')

  formated_runs = '\n'.join(considered_runs)
  print(f"These are the following runs gathered in the output:\n{formated_runs}\n")
  return considered_runs


def main():
  args = parser.parse_args()
  runs_path = args.runs_path
  output_path = args.output_path
  if not output_path:
    output_path = os.path.join(runs_path, 'all_pdbs')
  sequence_name = os.path.basename(os.path.realpath(runs_path))

  ignored_dir = ['all_pdbs', 'all_runs', 'msas', 'msas_colabfold']
  ignore_runs = [ run.replace('/', '') for run in args.runs_to_ignore ]
  ignored_dir.extend(ignore_runs)
  ignored_dir = [ directory for directory in sorted(list(set(ignored_dir)), key=str) if directory in os.listdir(runs_path) ]
  if ignored_dir:
    print(f"The following directories are ignored:\n{' - '.join(ignored_dir)}\n")

  if os.path.exists(output_path) and not args.only_ranking:
    print(f'{output_path} from previous iterations exists, deleting it then repeat.')
    rm(output_path)

  runs = check_all_runs(runs_path, ignored_dir)

  # create ranking json files
  delete_symlinks(runs_path, runs)
  create_symlink_without_ranked(runs_path, runs)
  rank_all(runs_path, runs, output_path) 
  if not args.only_ranking:
    pred_run_map = create_global_ranking(runs_path, runs, output_path) 
    print(f"Gathering {sequence_name}'s runs")
    if args.include_pickle:
      print("Pickle files are also included in the gathering.")
    move_and_rename(runs_path, pred_run_map, output_path, args.include_pickle)
  elif os.path.exists(output_path):
    rm(output_path)
  delete_symlinks(runs_path, runs)

if __name__ == "__main__":
  main()
