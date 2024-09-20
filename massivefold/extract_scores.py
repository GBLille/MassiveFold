#!/usr/bin/env python

import argparse
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import pickle
import json

parser = argparse.ArgumentParser(allow_abbrev=False)
parser.add_argument('--target_run', help='Path of the unfinished run you want to get the scores', required=True)
parser.add_argument('--verbose', choices=['ranking_confidence', 'ptm', 'iptm'], help='If you want scores to be displayed and which one')

def create_ranking(df_ranking, model_type='multimer'):
  ranking_metrics = [ ranking  for ranking in df_ranking.columns if ranking.startswith('ranking_') ]
  every_ranking = {}
  for metric in ranking_metrics:
    df = df_ranking.sort_values(metric, ascending=False).copy()
    score_type = metric.split('ranking_')[1]
    if score_type == 'debug' and model_type == 'multimer':
      score_type = 'iptm+ptm'
      df[score_type] = df[metric]
    else:
      df[score_type] = df[metric].apply(lambda x: float(x))

    scores = df.set_index('model')[score_type].to_dict()
    order = list(df['model'])
    ranking = {score_type: scores, 'order':order}
    every_ranking[metric] = ranking
  return every_ranking

def extract_scores(pkl_file):
  with open(pkl_file, 'rb') as pkl:
    data = pickle.load(pkl)
  keys = list(data.keys())
  
  if 'iptm' in keys: 
    return {'ranking_ptm': data['ptm'], 'ranking_iptm': data['iptm'], 'ranking_debug': data['ranking_confidence']}
  else:
    return {'ranking_ptm': data['ptm'], 'ranking_debug': data['ranking_confidence']}

def associate_pdb_pkl(path):
  pkl_presence = [ i for i in os.listdir(path) if i.endswith('.pkl') ]
  if pkl_presence:
    pkl_dir = path
  else:
    pkl_dir = os.path.join(path, './light_pkl')
  pdb = [ file for file in os.listdir(path) if file.endswith('.pdb') and f'result_model{file.split("model")[1].replace(".pdb", "")}.pkl' in os.listdir(pkl_dir) ]
  model = [ f'model{file.split("model")[1].replace(".pdb", "")}' for file in pdb  if f'result_model{file.split("model")[1].replace(".pdb", "")}.pkl' in os.listdir(pkl_dir) ]
  pkl = [ f'result_model{file.split("model")[1].replace(".pdb", "")}.pkl' for file in pdb  if f'result_model{file.split("model")[1].replace(".pdb", "")}.pkl' in os.listdir(pkl_dir) ]

  df = pd.DataFrame({'model': model, 'pdb': pdb, 'pkl': pkl})
  return df

def get_all_scores(df, path):
  pkl_presence = [ i for i in os.listdir(path) if i.endswith('.pkl') ]
  if pkl_presence:
    pkl_dir = path
  else:
    pkl_dir = os.path.join(path, './light_pkl')
  pkls = [ f"{pkl_dir}/{pkl}" for pkl in df['pkl'] ]
  results = list(map(extract_scores, pkls))
  """
  with ProcessPoolExecutor() as executor:
    results = list(executor.map(extract_scores, pkls))
  """
  keys = list(results[0].keys())
  df_scores = pd.DataFrame()
  for score_type in keys:
    scores = [ i[score_type] for i in results ]
    df_scores[score_type] = scores
  
  df_with_scores = pd.concat([df, df_scores], axis = 1)
  df_with_scores['model'] = df_with_scores['model'].apply(replace_version)
  return df_with_scores

def replace_version(row):
  if 'multimer' in row and 'v' not in row:
    parts = row.split('multimer')
    new_row = parts[0] + 'multimer_v1' + parts[1]
    return new_row
  return row

if __name__ == "__main__":
  args = parser.parse_args()
  df_association = associate_pdb_pkl(args.target_run)
  df_scored = get_all_scores(df_association, args.target_run)
  pdb = [ i for i in os.listdir(args.target_run) if i.endswith('.pdb') ]
  is_multimer = False
  for name in pdb:
    if 'multimer' in name:
      is_multimer = True    
  rankings = create_ranking(df_scored, is_multimer) 
  print(rankings)
  if args.verbose:
    score = f"ranking_{args.verbose}" if args.verbose != 'ranking_confidence' else 'ranking_debug'
    key = score.split('_')[1]
    key = key if key != 'debug' else 'iptm+ptm'
    print(json.dumps(rankings[score][key], indent=4))
  for ranking in rankings:
    with open(f'{args.target_run}/{ranking}.json', 'w') as json_ranking:
      json.dump(rankings[ranking], json_ranking, indent=4)
  
