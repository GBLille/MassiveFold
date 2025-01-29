#!/usr/bin/env python

import os
import json
import sys
import argparse
import numpy as np
import pandas as pd
from os import symlink
from shutil import copy  as cp, rmtree as rm

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
      True if "model" and "multimer" in pred
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
    for i, pred in enumerate(ordered_names):
      if "_unrelaxed_" in pred:
        full_filenames.append(f"ranked_{str(i)}_unrelaxed_{pred}.pdb")
      elif "_relaxed_" in pred:
        full_filenames.append(f"ranked_{str(i)}_relaxed_{pred}.pdb")
      else:
        do_not_exist.append(f"(ranked_{str(i)} => {pred})")
  elif is_alphafold3:
    full_filenames = [ f"ranked_{i}_{pred}.cif" for i, pred in enumerate(ordered_names) ]
  
  # check if these reconstituted files exist
  check_files_existence = lambda x: os.path.isfile(os.path.join(all_runs_path, run_name, x))
  do_not_exist.extend([ pred for pred in full_filenames if not check_files_existence(pred) ])
  assert not do_not_exist, f'Some files ({len(do_not_exist)}) for run {run_name} were not found: {", ".join(do_not_exist)}'

  return full_filenames

def rank_all(all_runs_path, all_runs, output_path, ranking_type="debug"):
  runs = [ os.path.join(all_runs_path, run) for run in all_runs ]
  ranking_key_score = ''
  all_models = pd.DataFrame()
  os.makedirs(output_path, exist_ok=True)

  for run in runs:
    single_run_models = pd.DataFrame()
    ranking_path = os.path.join(run, f'ranking_{ranking_type}.json')
    with open(ranking_path, 'r') as local_ranking_file:
      local_rank = json.load(local_ranking_file)
    if not ranking_key_score:
      ranking_key_score = list(local_rank.keys())[0]
    
    scores = list(local_rank[ranking_key_score].values())
    model_names = local_rank["order"]

    try:
      predictions = find_single_run_predictions(all_runs_path, run, model_names)
    except AssertionError as e:
      print(f"Assertion error: {str(e)[:300]}...")
      delete_symlinks(all_runs_path, all_runs)
      sys.exit()
    """
    reconstruct_prediction = lambda x, y: f"ranked_{y}_unrelaxed_{x}.pdb"
    predictions = list(map(reconstruct_prediction, model_names, range(len(model_names))))
    
    check_files_existence = lambda x: os.path.isfile(os.path.join(run, x))
  
    files_existence = list(map(check_files_existence, predictions))
    if not all(files_existence):
      indices = np.logical_not(files_existence).astype(int)
      not_found = list(np.array(predictions)[indices.astype(bool)])
      get_relaxed_version = lambda x: x.replace('unrelaxed', 'relaxed')
      relaxed_predictions = list(map(get_relaxed_version, not_found))
      relaxed_files_existence = list(map(check_files_existence, relaxed_predictions)) 
      if not all(relaxed_files_existence):
        secondary_indices = np.logical_not(relaxed_files_existence).astype(int)
        secondary_not_found = list(np.array(not_found)[secondary_indices.astype(bool)])
        print(f'/!\\ some predictions are not found in the run {run}, needs investigation')
        print(f"Predictions not present: {', '.join(secondary_not_found)}\n")
        delete_symlinks(all_runs_path, all_runs)
        sys.exit()
    """
    single_run_models['file'] = predictions
    single_run_models[ranking_key_score] = scores
    parameter_set = os.path.basename(os.path.normpath(run))
    single_run_models['parameters'] = parameter_set
    single_run_models["model_name"] = model_names
    #columns_order = ['parameters', 'file', ranking_key_score]
    #single_run_models = single_run_models[columns_order]
    all_models = pd.concat([all_models, single_run_models], axis=0) 
  
  all_models = all_models.sort_values(ranking_key_score, ascending=False, ignore_index=True)
  all_models['global_rank'] = all_models[ranking_key_score].rank(ascending=False, method='min').astype(int)
  columns_order = ['global_rank', ranking_key_score, 'parameters', 'file', "model_name"]
  all_models = all_models[columns_order]
  
  all_models_csv = all_models[ [col for col in all_models.columns if col != "model_name"] ]
  all_models_csv.to_csv(os.path.join(all_runs_path, 'ranking.csv'), index=None)
  print(all_models)
  return all_models

def move_and_rename(all_runs_path, run_names, output_path, ranking, do_include_pickles, do_include_rank):
  global_rank_order = ranking.to_dict(orient="records")
  models, runs, predictions, scores, mapped_names = [], [], [], [], []
  for i, prediction in enumerate(global_rank_order):
    prediction_run_name = prediction["parameters"]
    prediction_old_file = prediction["file"]
    model_name = prediction["model_name"]
    global_rank = prediction["global_rank"]
    ranking_score = prediction["ranking_score"]
    # copy the features if found
    if i == 0:
      features_old_name = os.path.join(all_runs_path, prediction_run_name, "features.pkl")
      features_new_name = os.path.join(output_path, "features.pkl")
      try:
        cp(features_old_name, features_new_name)
      except FileNotFoundError:
        print(f'features.pkl not found in {prediction_run_name} run')

    # copy the predictions
    old_pdb_path = os.path.join(all_runs_path, prediction_run_name, prediction_old_file)
    if old_pdb_path.endswith('.pdb'):
      structure_extension = 'pdb'
    elif old_pdb_path.endswith('.cif'):
      structure_extension = "cif"
    new_pdb_path = os.path.join(output_path, f"{model_name}.{structure_extension}")

    if do_include_rank:
      new_pdb_path = os.path.join(output_path, f"ranked_{global_rank - 1}_{prediction_run_name}_{model_name}.{structure_extension}")
      models.append(global_rank)
      runs.append(prediction_run_name)
      predictions.append(os.path.basename(old_pdb_path))
      scores.append(ranking_score)
      mapped_names.append(os.path.basename(new_pdb_path))
      
    cp(old_pdb_path, new_pdb_path)
    
    try:
      if do_include_pickles:
        pickles_path = os.path.join(all_runs_path, run_names[prediction])
        local_pred_name = f"model_{prediction.split('model_')[1]}" 
        if os.path.exists(os.path.join(pickles_path, 'light_pkl')):
          pickles_path = os.path.join(pickles_path, 'light_pkl')

        old_pdb_path = os.path.join(pickles_path, f"result_{local_pred_name}.pkl") 
        new_pdb_path = os.path.join(output_path, f"{prediction}.pkl")
        cp(old_pdb_path, new_pdb_path)
    except:
      print("--include_pickles should be fixed, not working anymore.")
  
  if do_include_rank:
    mapping = {"model": models, "run": runs, "prediction": predictions, "score": scores, "mapped_name": mapped_names}
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
  ranked_predictions = rank_all(runs_path, runs, output_path) 

  if not args.only_ranking:
    pred_run_map = create_global_ranking(runs_path, runs, output_path) 
    print(f"Gathering {sequence_name}'s runs")
    if args.include_pickles:
      print("Pickle files are also included in the gathering.")
    move_and_rename(runs_path, pred_run_map, output_path, ranked_predictions, args.include_pickles, args.include_rank)
  elif os.path.exists(output_path):
    rm(output_path)
  delete_symlinks(runs_path, runs)

if __name__ == "__main__":
  main()
