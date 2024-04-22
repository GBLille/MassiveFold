#!/usr/bin/env python

from absl import app, flags
import os
from os import symlink
import json
from shutil import copy  as cp, rmtree as rm
import sys

FLAGS = flags.FLAGS
flags.DEFINE_string(
  'runs_path',
  '',
  "Path containing the runs you want to gather.")
flags.DEFINE_list(
  'ignore_runs',
  [],
  "Coma separated list of runs to ignore in the results gathering.")

def create_global_ranking(all_runs_path, runs, ranking_type="debug"):
  map_pred_run = {}
  whole_ranking = {}
  runs = [ os.path.join(all_runs_path, run) for run in runs ]
  ranking_key_score = ''
  try:
    os.mkdir(f'{all_runs_path}/all_runs')
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
  with open(f"{all_runs_path}/all_runs/ranking_{ranking_type}.json", 'w') as fileout:
    fileout.write(json.dumps(global_ranking, indent=4)) 
  return map_pred_run

def move_and_rename(all_runs_path, run_names):
  with open(os.path.join(all_runs_path, 'all_runs', 'ranking_debug.json'), 'r') as rank_file:
    global_rank_order = json.load(rank_file)['order']
  for i, prediction in enumerate(global_rank_order):
    # copy the features
    if i == 0:
      features_old_name = os.path.join(all_runs_path, run_names[prediction], "features.pkl")
      features_new_name = os.path.join(all_runs_path, 'all_runs', "features.pkl")
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
    new_pdb_path = os.path.join(all_runs_path, 'all_runs', f"ranked_{i}_{prediction}.pdb")
    cp(old_pdb_path, new_pdb_path)

def create_symlink_without_ranked(all_runs_path, runs):
  runs = [ os.path.join(all_runs_path, run) for run in runs ]
  for run in runs:
    path = os.path.join(all_runs_path, run)
    ranked_preds = [pred for pred in os.listdir(path) if pred.startswith('ranked_')]
    for pred in ranked_preds:
      old_name = os.path.join(path, pred)
      new_name = os.path.join(path, pred.split('_', 2)[-1])
      try:
        os.symlink(old_name, new_name)
      except FileExistsError:
        pass

def delete_symlinks(all_runs_path, runs):
  runs = [ os.path.join(all_runs_path, run) for run in runs ]
  for run in runs:
    path_run = os.path.join(all_runs_path, run)
    files = os.listdir(path_run)
    for symlink in files:
      symlink_file = os.path.join(path_run, symlink)
      if os.path.islink(symlink_file):
        os.remove(symlink_file)

def check_all_runs(all_runs_path, ignored_runs, ignored_directories, ranking_type='debug'):
  considered_runs = []
  
  for run in os.listdir(all_runs_path):
    if run not in ignored_directories and run not in ignored_runs:
      ranking_path = os.path.join(all_runs_path, run, f'ranking_{ranking_type}.json')
      if os.path.isfile(ranking_path):
        try:
          with open(ranking_path, 'r') as ranking_file:
            json.load(ranking_file)
          considered_runs.append(run)
        except:
          print(f'Something went wrong with ranking_debug.json for run {run}')

  formated_runs = '\n'.join(considered_runs)
  print(f"These are the following runs gathered in all_runs:\n{formated_runs}")
  return considered_runs

def main(argv):
  sequence_name = os.path.basename(FLAGS.runs_path)
  runs_path = os.path.abspath(FLAGS.runs_path)
  ignored_runs = FLAGS.ignore_runs
  ignored_directories = ['all_runs', 'msas', 'msas_colabfold']
 
  all_runs = os.path.join(runs_path, 'all_runs')
  if os.path.exists(all_runs):
    print('all_runs from previous iterations exists, deleting it then repeat.')
    rm(all_runs)
  runs = check_all_runs(runs_path, ignored_runs, ignored_directories)
  # create ranking json files
  delete_symlinks(runs_path, runs)
  create_symlink_without_ranked(runs_path, runs)
  pred_run_map = create_global_ranking(runs_path, runs) 
  move_and_rename(runs_path, pred_run_map)
  delete_symlinks(runs_path, runs)

if __name__ == "__main__":
  app.run(main)
