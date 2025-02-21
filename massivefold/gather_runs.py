#!/usr/bin/env python

import os
import json
import sys
import argparse
import numpy as np
import pandas as pd
from os import symlink
from shutil import copy as cp, rmtree as rm

parser = argparse.ArgumentParser(allow_abbrev=False)
parser.add_argument('--runs_path', help='Path to the runs you want to gather', required=True)
parser.add_argument('--ignore', dest='runs_to_ignore', nargs='+', help="List of run names or" 
" subdirectories in the run path set in --runs_path. All these runs (separated by a single space)" 
" won't be taken into account in the runs gathering", default=[])
parser.add_argument('--output_path', help="Path to the output of the runs gathering (default: <RUNS_PATH>/all_pdbs)", required=False)
parser.add_argument('--only_ranking', help="Skips the run gathering, only output a csv ranking file", required=False, action='store_true')
parser.add_argument('--include_pickles', help="If specified, include the .pkl pickle" 
" files in the gathered results stored in <OUTPUT_PATH>", required=False, action='store_true')
parser.add_argument('--include_rank', help="Include a ranking (ambiguous for ties resolving)" 
"in the name of the ranked files, also provides a mapping from original files to the gathered files.", required=False, action='store_true')

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

def find_single_run_predictions(all_runs_path: str, run_name: str, ordered_names: list):
  """
  Given a MassiveFold run directory and the names of this run (ordered by score but
  without "ranked_{n}_, find the full filenames. In case some of the predictions are
  not found, interrupt the script while cleaning the symlink that may have been created.
  """
  # find run's nature (AF2/AF3)
  run_name = os.path.basename(run_name)
  is_alphafold2 = False
  is_alphafold3 = False
  names_contain_NN_model_names = [
      True if "model" and ("multimer" in pred or "ptm" in pred)
      else False
      for pred in ordered_names
 ]
  if all(names_contain_NN_model_names):
    is_alphafold2 = True
  names_contain_seed = [True if "seed" and "sample" in pred else False for pred in ordered_names]
  if all(names_contain_seed):
    is_alphafold3 = True
  assert not (is_alphafold2 and is_alphafold3), \
  "Could not determinate if the run was AlphaFold2 or AlphaFold3 based, detected as both."
  assert is_alphafold2 or is_alphafold3, \
  "Could not determinate if the run was AlphaFold2 or AlphaFold3 based, detected as neither."

  # reconstitute full filenames
  do_not_exist = []
  full_filenames = []
  if is_alphafold2:
    all_preds = os.listdir(os.path.join(all_runs_path, run_name))
    for i, pred in enumerate(ordered_names):
      prefix = f"ranked_{str(i)}_"
      suffix = f"{pred}.pdb"
      matches = [ i for i in all_preds if i.startswith(prefix) and i.endswith(suffix) ]
      if len(matches) == 0 or len(matches) > 1:
        do_not_exist.append(f"(ranked_{str(i)} => {pred})")
      else:
        full_filenames.append(matches[0])
  elif is_alphafold3:
    full_filenames = [ f"ranked_{i}_{pred}.cif" for i, pred in enumerate(ordered_names) ]
  
  # check if these reconstituted files exist
  does_file_exist = lambda x: os.path.isfile(os.path.join(all_runs_path, run_name, x))
  do_not_exist.extend([ pred for pred in full_filenames if not does_file_exist(pred) ])
  assert not do_not_exist, f'Some files ({len(do_not_exist)}) for run {run_name} were not found: {", ".join(do_not_exist)}'

  return full_filenames

def rank_all(all_runs_path, all_runs, output_path, ranking_type="debug"):
  runs = [ os.path.join(all_runs_path, run) for run in all_runs ]

  all_models = pd.DataFrame()
  os.makedirs(output_path, exist_ok=True)

  for run in runs:
    single_run_models = pd.DataFrame()
    ranking_path = os.path.join(run, f'ranking_debug.json')
    ranking_score_path = os.path.join(run, f'ranking_{ranking_type}.json')
    with open(ranking_path, 'r') as local_ranking_file:
      local_rank = json.load(local_ranking_file)

    model_names = local_rank["order"]
    score_ranking = json.load(open(ranking_score_path, 'r'))
    ranking_key_score = list(score_ranking.keys())[0]
    scores = [ score_ranking[ranking_key_score][model] for model in model_names ]
    try:
      predictions = find_single_run_predictions(all_runs_path, run, model_names)
    except AssertionError as e:
      print(f"Assertion error: {str(e)[:300]}...")
      delete_symlinks(all_runs_path, all_runs)
      sys.exit()
    single_run_models['file'] = predictions
    single_run_models[ranking_key_score] = scores
    parameter_set = os.path.basename(os.path.normpath(run))
    single_run_models['parameters'] = parameter_set
    single_run_models["model_name"] = model_names
    all_models = pd.concat([all_models, single_run_models], axis=0) 
  
  all_models = all_models.sort_values(ranking_key_score, ascending=False, ignore_index=True)
  all_models['global_rank'] = all_models[ranking_key_score].rank(ascending=False, method='min').astype(int)
  columns_order = ['global_rank', ranking_key_score, 'parameters', 'file', "model_name"]
  all_models = all_models[columns_order]
  
  all_models_csv = all_models[ [col for col in all_models.columns if col != "model_name"] ]
  return all_models

def move_and_rename(all_runs_path, run_names, output_path, ranking, do_include_pickles, do_include_rank):
  score_key = ranking.columns[1]
  global_rank_order = ranking.to_dict(orient="records")
  models, runs, predictions, scores, mapped_names = [], [], [], [], []
  for i, prediction in enumerate(global_rank_order):
    run_name = prediction["parameters"]
    prediction_old_file = prediction["file"]
    model_name = prediction["model_name"]
    global_rank = prediction["global_rank"]
    ranking_score = prediction[score_key]
    # copy the features if found
    if i == 0:
      features_old_name = os.path.join(all_runs_path, run_name, "features.pkl")
      features_new_name = os.path.join(output_path, "features.pkl")
      try:
        cp(features_old_name, features_new_name)
      except FileNotFoundError:
        print(f'features.pkl not found in {run_name} run')

    # copy the predictions
    old_pdb_path = os.path.join(all_runs_path, run_name, prediction_old_file)
    if old_pdb_path.endswith('.pdb'):
      structure_extension = 'pdb'
    elif old_pdb_path.endswith('.cif'):
      structure_extension = "cif"
    pdb_file = f"{run_name}_{model_name}.{structure_extension}"

    if do_include_rank:
      pdb_file = f"ranked_{global_rank - 1}_" + pdb_file
      models.append(global_rank)
      runs.append(run_name)
      predictions.append(os.path.basename(old_pdb_path))
      scores.append(ranking_score)
      mapped_names.append(os.path.basename(pdb_file))
      
    new_pdb_path = os.path.join(output_path, pdb_file)
    cp(old_pdb_path, new_pdb_path)
    
    if do_include_pickles:
      prediction_name = prediction["parameters"]  + "_" + prediction["model_name"]
      pickles_path = os.path.join(all_runs_path, run_names[prediction_name])
      if os.path.exists(os.path.join(pickles_path, 'light_pkl')):
        pickles_path = os.path.join(pickles_path, 'light_pkl')

      old_pdb_path = os.path.join(pickles_path, f"result_{prediction['model_name']}.pkl") 
      new_pdb_path = os.path.join(output_path, f"{prediction_name}.pkl")
      cp(old_pdb_path, new_pdb_path)
  
  if do_include_rank:
    mapping = {"model": models, "run": runs, "prediction": predictions, score_key: scores, "mapped_name": mapped_names}
    pd.DataFrame(mapping).to_csv(os.path.join(output_path, 'map.csv'), index=False)

def create_symlink_without_ranked(all_runs_path, runs):
  for run in runs:
    path = os.path.join(all_runs_path, run)
    ranked_preds = [ pred for pred in os.listdir(path) if pred.startswith('ranked_') ]
    # to remove when unrelaxed is implemented
    corrected_preds = [ pred.replace('relaxed', 'unrelaxed') if '_relaxed' in pred else pred for pred in ranked_preds ]
    for pred, corr in zip(ranked_preds, corrected_preds):
      old_name = os.path.realpath(os.path.join(path, pred))
      new_name = os.path.join(path, corr.split('_', 2)[-1])
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

  ranking_types = ["debug", "iptm", "ptm", "plddt"]
  ignored_dir = ['all_pdbs', 'all_runs', 'msas', 'msas_colabfold', "msas_alphafold3"]
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

  for i, ranking_type in enumerate(ranking_types):
    try:
      ranked_predictions = rank_all(runs_path, runs, output_path, ranking_type) 
    except FileNotFoundError as e:
      print(f"No ranking for {ranking_type} metric.")
    if i == 0:
      whole_prediction_ranking = ranked_predictions
    else:
      new_metric = [ column_new for column_new in ranked_predictions.columns if column_new not in whole_prediction_ranking.columns ]
      whole_prediction_ranking[new_metric] = ranked_predictions[new_metric]

  assert not whole_prediction_ranking.empty
  print(whole_prediction_ranking)
  whole_prediction_ranking.to_csv(os.path.join(runs_path, 'ranking.csv'), index=None)
  whole_prediction_ranking.to_csv(os.path.join(output_path, 'ranking.csv'), index=None)

  if not args.only_ranking:
    pred_run_map = create_global_ranking(runs_path, runs, output_path) 
    print(f"Gathering {sequence_name}'s runs")
    if args.include_pickles:
      print("Pickle files are also included in the gathering.")
    move_and_rename(runs_path, pred_run_map, output_path, whole_prediction_ranking, args.include_pickles, args.include_rank)
  elif os.path.exists(output_path):
    rm(output_path)
  delete_symlinks(runs_path, runs)

if __name__ == "__main__":
  main()
