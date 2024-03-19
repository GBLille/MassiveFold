#!/usr/bin/env python

from absl import app, flags
import os
from os import symlink
import json
from shutil import copy  as cp
import sys

FLAGS = flags.FLAGS
flags.DEFINE_string(
    'runs_path',
    '',
    "Path containing the runs you want to gather.")

def create_global_ranking(all_runs_path, ranking_type="debug"):
  map_pred_run = {}
  whole_ranking = {}
  ranking_key_score = ''

  try:
    os.mkdir(f'{all_runs_path}/all_runs')
  except FileExistsError:
    pass
  for run in os.listdir(FLAGS.runs_path):
    if run != 'all_runs' and run != 'msas':
      ranking_path = os.path.join(all_runs_path, run, f'ranking_{ranking_type}.json')
      with open(ranking_path, 'r') as local_ranking_file:
        local_rank = json.load(local_ranking_file)
        if not ranking_key_score:
          ranking_key_score = list(local_rank.keys())[0]
      local_rank[ranking_key_score] = \
      {f"{run}_{pred}": local_rank[ranking_key_score][pred] for pred in local_rank[ranking_key_score]}
      local_rank['order'] = [f"{run}_{pred}" for pred in local_rank['order']]
      map_pred_run.update({pred: run for pred in local_rank['order']})
      whole_ranking.update(local_rank[ranking_key_score])
  sorted_predictions = sorted(whole_ranking.items(), key=lambda x:x[1], reverse=True)
  iptm_ptm = dict(sorted_predictions)
  order = sorted(whole_ranking, reverse=True, key=whole_ranking.get)
  global_ranking = {ranking_key_score: iptm_ptm, 'order': order}
  with open(f"{FLAGS.runs_path}/all_runs/ranking_{ranking_type}.json", 'w') as fileout:
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
 
def create_symlink_without_ranked(all_runs_path):
  runs = os.listdir(all_runs_path)
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

def delete_symlinks(all_runs_path):
  runs = os.listdir(all_runs_path)
  for run in runs:
    path_run = os.path.join(all_runs_path, run)
    files = os.listdir(path_run)
    for symlink in files:
      symlink_file = os.path.join(path_run, symlink)
      if os.path.islink(symlink_file):
        os.remove(symlink_file)

def main(argv):
  FLAGS.runs_path = os.path.abspath(FLAGS.runs_path)
  sequence_name = os.path.basename(FLAGS.runs_path)
  
  # create ranking json files
  delete_symlinks(FLAGS.runs_path)
  create_symlink_without_ranked(FLAGS.runs_path)
  pred_run_map = create_global_ranking(FLAGS.runs_path) 
  move_and_rename(FLAGS.runs_path, pred_run_map)
  delete_symlinks(FLAGS.runs_path)

if __name__ == "__main__":
  app.run(main)     
