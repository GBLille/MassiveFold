#!/usr/bin/env python

import argparse
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import pickle
import json

parser = argparse.ArgumentParser(allow_abbrev=False)
parser.add_argument('--target_run', help='Path of the unfinished run you want to get the scores', required=True)


def create_ranking(df_ranking):
  ranking_metrics = [ ranking  for ranking in df_ranking.columns if ranking.startswith('ranking_') ]
  every_ranking = {}
  for metric in ranking_metrics:
    df = df_ranking.sort_values(metric, ascending=False).copy()
    score_type = metric.split('ranking_')[1]
    if score_type == 'debug':
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
  
  pdb = [ file for file in os.listdir(path) if file.startswith('unrelaxed_') if f'result_{file.split("unrelaxed_")[1].split(".pdb")[0]}.pkl' in os.listdir(path) ]
  model = [ f'{file.split("unrelaxed_")[1].split(".pdb")[0]}' for file in pdb  if f'result_{file.split("unrelaxed_")[1].split(".pdb")[0]}.pkl' in os.listdir(path) ]
  pkl = [ f'result_{file.split("unrelaxed_")[1].split(".pdb")[0]}.pkl' for file in pdb  if f'result_{file.split("unrelaxed_")[1].split(".pdb")[0]}.pkl' in os.listdir(path) ]

  df = pd.DataFrame({'model': model, 'pdb': pdb, 'pkl': pkl})
  return df

def get_all_scores(df, path):
  pkls = [ f"{path}/{pkl}" for pkl in df['pkl'] ]
  with ProcessPoolExecutor() as executor:
    results = list(executor.map(extract_scores, pkls))
  keys = list(results[0].keys())
  df_scores = pd.DataFrame()
  for score_type in keys:
    scores = [ i[score_type] for i in results ]
    df_scores[score_type] = scores
  
  df_with_scores = pd.concat([df, df_scores], axis = 1)
  return df_with_scores

if __name__ == "__main__":
  args = parser.parse_args()
  
  df_association = associate_pdb_pkl(args.target_run)
  df_scored = get_all_scores(df_association, args.target_run)
  rankings = create_ranking(df_scored) 
  
  for ranking in rankings:
    with open(f'{args.target_run}/{ranking}.json', 'w') as json_ranking:
      json.dump(rankings[ranking], json_ranking, indent=4)
  
