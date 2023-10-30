#!/usr/bin/env python

from absl import app, flags
import os
import json
from shutil import copy as cp, rmtree as rm, move as mv

FLAGS = flags.FLAGS
flags.DEFINE_string('batches_path', '', "Path of all batches containing the ranking files to add in the global ranking.")
flags.DEFINE_string('global_ranking_path', 'global_ranking.json', "Global ranking file name, in the 'global path' directory.")

def create_global_ranking(all_batches_path, jobname, return_map=False):
  map_pred_batch = {}
  whole_ranking = {}
  for batch in os.listdir(FLAGS.batches_path):
    if batch.startswith('batch'):
      ranking_path = os.path.join(FLAGS.batches_path, batch, jobname, 'ranking_debug.json')
      with open(ranking_path, 'r') as local_ranking_file:
        local_rank = json.load(local_ranking_file)
      map_pred_batch.update({pred: batch for pred in local_rank['order']})
      whole_ranking.update(local_rank["iptm+ptm"])
  sorted_predictions = sorted(whole_ranking.items(), key=lambda x:x[1], reverse=True)
  iptm_ptm = dict(sorted_predictions)
  order = sorted(whole_ranking, reverse=True, key=whole_ranking.get)
  global_ranking = {'iptm+ptm': iptm_ptm, 'order': order}
  with open(f"{FLAGS.batches_path}/ranking_debug.json", 'w') as fileout:
    fileout.write(json.dumps(global_ranking, indent=4)) 
  return map_pred_batch

def move_and_rename(all_batches_path, pred_batch_map, jobname):
  with open(os.path.join(all_batches_path, 'ranking_debug.json'), 'r') as rank_file:
    global_rank_order = json.load(rank_file)['order']
  for i, prediction in enumerate(global_rank_order):
    # copy the features
    if i == 0:
      cp(os.path.join(all_batches_path, pred_batch_map[prediction], jobname, "features.pkl"), os.path.join(all_batches_path, "features.pkl"))
    
    # copy the predictions
    if os.path.exists(os.path.join(all_batches_path, pred_batch_map[prediction], jobname, f"relaxed_{prediction}.pdb")):
      pred_new_name = f"relaxed_{prediction}.pdb"
    else:
      pred_new_name = f"unrelaxed_{prediction}.pdb"
    
    # Move pdb files and rename with rank
    old_pdb_path = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, pred_new_name) 
    new_pdb_path = os.path.join(all_batches_path, f"ranked_{i}_{pred_new_name}")
    cp(old_pdb_path, new_pdb_path)
    
    # Move pkl files
    pkl_name = f"result_{prediction}.pkl"
    old_pkl_path = os.path.join(all_batches_path, pred_batch_map[prediction], jobname, pkl_name)
    new_pkl_path = os.path.join(all_batches_path, pkl_name)
    mv(old_pkl_path, new_pkl_path)

def remove_batch_dirs(all_batches_path):
  batch_dirs = [d for d in os.listdir(all_batches_path) if d.startswith('batch')]
  for batch_dir in batch_dirs:
    rm(os.path.join(all_batches_path, batch_dir))
    
def main(argv):
  FLAGS.batches_path = os.path.abspath(FLAGS.batches_path)
  sequence_name = os.path.basename(FLAGS.batches_path)

  # check if batches parent directory is run_name
  if os.path.dirname(FLAGS.batches_path) != "output_array":
    sequence_name = os.path.basename(os.path.dirname(FLAGS.batches_path))

  pred_batch_map = create_global_ranking(FLAGS.batches_path, sequence_name)
  move_and_rename(FLAGS.batches_path, pred_batch_map, sequence_name)
  remove_batch_dirs(FLAGS.batches_path)

if __name__ == "__main__":
  app.run(main)     
