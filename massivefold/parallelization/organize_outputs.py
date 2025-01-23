#!/usr/bin/env python

from absl import app, flags
import os
import json
from shutil import copy as cp, rmtree as rm, move as mv
import sys

FLAGS = flags.FLAGS
flags.DEFINE_string('batches_path', '', "Path of all batches containing the ranking files to add in the global ranking.")

def create_global_ranking(all_batches_path, jobname, ranking_type="debug"):
  map_pred_batch = {}
  whole_ranking = {}
  ranking_key_score = ''
  for batch in os.listdir(all_batches_path):
    if batch.startswith('batch'):
      ranking_path = os.path.join(all_batches_path, batch, jobname, f'ranking_{ranking_type}.json')
      with open(ranking_path, 'r') as local_ranking_file:
        local_rank = json.load(local_ranking_file)
        if not ranking_key_score:
          ranking_key_score = list(local_rank.keys())[0]
      map_pred_batch.update({pred: batch for pred in local_rank['order']})
      whole_ranking.update(local_rank[ranking_key_score])
  sorted_predictions = sorted(whole_ranking.items(), key=lambda x:x[1], reverse=True)
  iptm_ptm = dict(sorted_predictions)
  order = sorted(whole_ranking, reverse=True, key=whole_ranking.get)
  global_ranking = {ranking_key_score: iptm_ptm, 'order': order}
  with open(f"{all_batches_path}/ranking_{ranking_type}.json", 'w') as fileout:
    fileout.write(json.dumps(global_ranking, indent=4)) 

  return map_pred_batch

def move_and_rename(all_batches_path, pred_batch_map, jobname):
  with open(os.path.join(all_batches_path, 'ranking_debug.json'), 'r') as rank_file:
    global_rank_order = json.load(rank_file)['order']
  for i, prediction in enumerate(global_rank_order):
    # copy the features
    if i == 0:
      features = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, "features.pkl")
      if os.path.isfile(features):
        cp(features, os.path.join(all_batches_path, "features.pkl"))
      else:
        print('Either using colabfold or error encountered while copying features.pkl')
      files = os.path.join(all_batches_path, pred_batch_map[prediction], jobname)
      for file in os.listdir(files):
        if file.endswith('coverage.png'):
          coverage_plot = os.path.join(files, file)
          os.mkdir(os.path.join(all_batches_path, "./plots"))
          cp(coverage_plot, os.path.join(all_batches_path, "./plots/alignment_coverage.png"))
      else:
        print('Either not using colabfold, or coverage plot not found')

    # copy the predictions
    if os.path.exists(os.path.join(all_batches_path, pred_batch_map[prediction], jobname, f"relaxed_{prediction}.pdb")):
      pred_new_name = f"relaxed_{prediction}.pdb"
    elif os.path.exists(os.path.join(all_batches_path, pred_batch_map[prediction], jobname, f"unrelaxed_{prediction}.pdb")):
      pred_new_name = f"unrelaxed_{prediction}.pdb"
    else:
      pred_new_name = f"{prediction}.cif"
    
    # Move confidence files (chain iptm etc)
    confidence_path = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, 'confidences')
    if os.path.exists(confidence_path):
      old_confidence_path = os.path.join(
        all_batches_path,
        pred_batch_map[prediction],
        jobname,
        'confidences',
        f"{pred_new_name.replace('.cif', '.json')}") 
      global_confidence_path = os.path.join(all_batches_path, 'confidences')
      if not os.path.exists(global_confidence_path):
        os.makedirs(global_confidence_path)
      new_confidence_path = os.path.join(global_confidence_path, f"ranked_{i}_{pred_new_name.replace('.cif', '.json')}")
      try:
        mv(old_confidence_path, new_confidence_path)
      except FileNotFoundError as e:
        print(e)
        print(f"{pred_batch_map[prediction]}/confidences/{pred_new_name.replace('.cif', '.json')} not found")

    # Move pdb files and rename with rank
    old_pdb_path = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, pred_new_name) 
    new_pdb_path = os.path.join(all_batches_path, f"ranked_{i}_{pred_new_name}")
    try:
      mv(old_pdb_path, new_pdb_path)
    except FileNotFoundError as e:
      print(e)
      print(f"{pred_batch_map[prediction]}/ranked_{i}_{pred_new_name} does not exist, probably score < --min_score.")

    # Move pkl files
    pkl_name = f"result_{prediction}.pkl"
    old_pkl_path = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, pkl_name)
    new_pkl_path = os.path.join(all_batches_path, pkl_name)
    try:
      mv(old_pkl_path, new_pkl_path)
    except FileNotFoundError:
      print(f"{pred_batch_map[prediction]}/result_{prediction}.pkl does not exist, probably score < --min_score.")

def remove_batch_dirs(all_batches_path):
  batch_dirs = [d for d in os.listdir(all_batches_path) if d.startswith('batch')]
  for batch_dir in batch_dirs:
    rm(os.path.join(all_batches_path, batch_dir))

def main(argv):
  batches_path = os.path.normpath(FLAGS.batches_path)
  sequence_name = os.path.basename(os.path.dirname(batches_path))
  run_name = os.path.basename(batches_path)
  batch_0 = os.path.join(batches_path, "batch_0")

  if os.path.exists(os.path.join(batch_0, "ranking_debug.json")):
    print("AlphaFold3 output detected")
    try:
      mv(os.path.join(batch_0, "TERMS_OF_USE.md"), batches_path)
    except:
      pass
    all_batches = [ i for i in os.listdir(batches_path) if i.startswith('batch_') and os.path.isdir(os.path.join(batches_path, i)) ]
    for batch in all_batches:
      batch_dir = os.path.join(batches_path, batch)
      new_location = os.path.join(batch_dir, sequence_name)
      os.makedirs(new_location)
      output_files = [ os.path.join(batch_dir, i) for i in os.listdir(batch_dir) if i.endswith('.json') or i.endswith('.cif') or i.endswith('.pkl') ]

      confidence_path = os.path.join(batch_dir, 'confidences')
      new_confidence_path = os.path.join(new_location, 'confidences')

      os.makedirs(new_confidence_path)
      for confidence_file in os.listdir(confidence_path):
        cp(os.path.join(confidence_path, confidence_file), new_confidence_path)
      for file in output_files: 
        new_file = os.path.join(new_location, os.path.basename(file))
        mv(file, new_file)

  # create ranking json files
  pred_batch_map = create_global_ranking(batches_path, sequence_name)
  if os.path.isfile(f"{batches_path}/batch_0/{sequence_name}/ranking_iptm.json"):
    create_global_ranking(batches_path, sequence_name, 'iptm')
  if os.path.isfile(f"{batches_path}/batch_0/{sequence_name}/ranking_ptm.json"):
    create_global_ranking(batches_path, sequence_name, 'ptm')

  # organize output directory
  move_and_rename(batches_path, pred_batch_map, sequence_name)
  remove_batch_dirs(batches_path)

if __name__ == "__main__":
  app.run(main)     
