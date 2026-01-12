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
parser.add_argument('--include_rank', help="Include a ranking (ambiguous for ties resolving) on 0.8 x iptm + 0.2 x ptm" 
"in the name of the ranked files, also provides a mapping from original files to the gathered files.", required=False, action='store_true')

def global_rank_to_json(ranking, output_path):
  map_pred_run = dict(zip(list(ranking["parameters"] + "_" + ranking["model_name"]), list(ranking['parameters'])))
  info = ['global_rank', 'parameters', 'file', "model_name"]
  for metric in ranking.columns:
    if metric in info:
      continue

    ranking_type = metric
    if metric == 'iptm+ptm' or metric == 'ranking_score':
      ranking_type = "debug"
    df = ranking[~ranking[metric].isna()].copy()
    df = df.sort_values(metric, ascending=False, ignore_index=True)
    order = list(df["parameters"] + "_" + df["model_name"])
    model_scores = dict(zip(order, list(df[metric])))
    global_ranking = {metric: model_scores, "order": order}
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

  ranked_per_run = {}
  for run in runs:
    single_run_models = pd.DataFrame()
    ranking_path = os.path.join(run, f'ranking_debug.json')
    run_score_path = os.path.join(run, f'ranking_{ranking_type}.json')
    with open(ranking_path, 'r') as local_ranking_file:
      local_rank = json.load(local_ranking_file)

    model_names = local_rank["order"]
    score_ranking = json.load(open(run_score_path, 'r'))
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
    ranked_per_run[os.path.basename(run)] = single_run_models
    #all_models = pd.concat([all_models, single_run_models], axis=0)
  
  return ranked_per_run

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

      old_pkl_path = os.path.join(pickles_path, f"result_{prediction['model_name']}.pkl")
      new_pkl_path = os.path.join(output_path, f"{prediction_name}.pkl")
      try:
        cp(old_pkl_path, new_pkl_path)
      except:
        pass

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

  formated_runs = ' - '.join(considered_runs)
  print(f"These are the {len(considered_runs)} following runs gathered in the output:\n{formated_runs}\n")
  return considered_runs

def create_global_ranking(runs, runs_path, output_path, ranking_types):
  all_metrics_ranking = []
  for i, ranking_type in enumerate(ranking_types):
    try:
      ranking_per_run = rank_all(runs_path, runs, output_path, ranking_type) 
      all_metrics_ranking.append(ranking_per_run)
    except FileNotFoundError as e:
      print(f"No ranking for {ranking_type} metric.")
  all_run_names = [ list(metric.keys()) for metric in all_metrics_ranking ]
  run_names_homogeneity = [ run_names == list(all_metrics_ranking[0].keys()) for run_names in all_run_names ]
  assert all(run_names_homogeneity), \
  f"Not the same runs found: {', '.join([str(i) for i in all_run_names])}"
  run_names = list(all_metrics_ranking[0].keys())

  columns_to_merge_on = ['file', 'parameters', 'model_name']
  score_keys = []
  all_runs_ranking = pd.DataFrame()
  for run in run_names:
    run_ranking = all_metrics_ranking[0][run]
    run_key_score = [ col for col in run_ranking.columns if col not in columns_to_merge_on ][0]
    for i in range(1, len(all_metrics_ranking)):
      run_ranking = run_ranking.merge(all_metrics_ranking[i][run], on=columns_to_merge_on, how='left')
    if run_key_score == "ranking_score":
      run_key_score = "af3_ranking_score"

      if "iptm" in run_ranking.columns:
        run_ranking["iptm+ptm"] = 0.8*run_ranking["iptm"] + 0.2*run_ranking["ptm"]

      run_ranking = run_ranking.rename(columns={"ranking_score": "af3_ranking_score"})
    all_runs_ranking = pd.concat([all_runs_ranking, run_ranking], axis=0)
    score_keys.append(run_key_score)

  if len(set(score_keys)) == 1:
    common_key_score = list(set(score_keys))[0]
    ordering_score = common_key_score
  elif "af3_ranking_score" in score_keys and "iptm+ptm" in score_keys:
    common_key_score = "iptm+ptm"
    ordering_score = [common_key_score, 'af3_ranking_score']
  elif "af3_ranking_score" in score_keys and "plddts" in score_keys:
    #common_key_score = "plddts"
    common_key_score = "ptm" #temporary fix before ranking on plddts is implemented for af3
    ordering_score = [common_key_score, 'af3_ranking_score']
  else:
    raise ValueError(f"ranking_debug.json in some runs uses different metrics: {score_keys}")

  all_runs_ranking = all_runs_ranking.sort_values(ordering_score, ascending=False, ignore_index=True)
  all_runs_ranking["global_rank"] = all_runs_ranking[common_key_score].rank(ascending=False, method='min').astype(int)
  columns_order = ['global_rank', common_key_score, 'parameters', 'file', "model_name"]
  columns_order.extend([col for col in all_runs_ranking.columns if col not in columns_order])
  whole_prediction_ranking = all_runs_ranking[columns_order]
  return whole_prediction_ranking

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
  assert runs, "There should be at least one run to gather"
  delete_symlinks(runs_path, runs)
  create_symlink_without_ranked(runs_path, runs)
  whole_prediction_ranking = create_global_ranking(runs, runs_path, output_path, ranking_types)

  assert not whole_prediction_ranking.empty
  other_csv = [ csv for csv in os.listdir(runs_path) if csv.endswith('.csv') and csv != "ranking.csv"]
  for csv in other_csv:
    csv_file = os.path.join(runs_path, csv)
    df = pd.read_csv(csv_file)
    if "i-plddt" in df.columns:
      print(f"DataFrame containing i-plddt values: {csv_file}")
      print(df)
      whole_prediction_ranking["extension"] = whole_prediction_ranking["file"].str.split('.').str[1]
      whole_prediction_ranking["Models"] = whole_prediction_ranking["parameters"] + '_' + whole_prediction_ranking["model_name"] + '.' + whole_prediction_ranking["extension"]
      whole_prediction_ranking = pd.merge(whole_prediction_ranking, df, on="Models", how='left')
      whole_prediction_ranking = whole_prediction_ranking.drop(columns={"Models", "extension"})
      break

  print(whole_prediction_ranking)
  whole_prediction_ranking.to_csv(os.path.join(runs_path, 'ranking.csv'), index=None)
  whole_prediction_ranking.to_csv(os.path.join(output_path, 'ranking.csv'), index=None)

  if not args.only_ranking:
    # create ranking json files
    pred_run_map = global_rank_to_json(whole_prediction_ranking, output_path)
    print(f"Gathering {sequence_name}'s runs")
    if args.include_pickles:
      print("Pickle files are also included in the gathering.")
    move_and_rename(runs_path, pred_run_map, output_path, whole_prediction_ranking, args.include_pickles, args.include_rank)
  elif os.path.exists(output_path):
    rm(output_path)
  delete_symlinks(runs_path, runs)

if __name__ == "__main__":
  main()
