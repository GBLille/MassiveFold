#!/usr/bin/env python
import os
from absl import flags, app
from Bio import SeqIO
import pickle
from shutil import move as mv
import pandas as pd
import json
import sys
import shutil

FLAGS = flags.FLAGS
flags.DEFINE_enum(
  'conversion',
  None,
  ['input', 'output'],
  "What to convert."
  "'input' to get the fasta format of colabfold from a traditionnal multichain pdb."
  "'output' to transform colabfold output to the alphafold one.")
flags.DEFINE_string(
  'to_convert',
  '',
  'Path of the fasta file (for --conversion=input) or output directory (for --conversion=output)'
  'To convert')
flags.DEFINE_string(
  'batches_file',
  '',
  'Path to batches file. If --conversion=output, this file is necessary.')
flags.DEFINE_bool(
  'do_rename',
  True,
  'To rename file or not')


def convert_fasta(fasta_path:str):
  records = list(SeqIO.parse(fasta_path, "fasta"))
  fasta_dir = os.path.dirname(fasta_path)
  fasta_file = os.path.basename(fasta_path).split('.')[0]

  converted = ''
  ids = '>'
  sequences = ''
  for record in records:
    ids += f"{record.id}:"
    sequences += f"{record.seq}:"

  converted = f"{ids[:-1]}\n{sequences[:-1]}"
  
  with open(f"{fasta_dir}/converted_for_colabfold/{fasta_file}.fasta", 'w') as output:
    output.write(converted)

def rename_pkl(pkl_files:list, output_path:str, pred_shift:int, sep:str):
  extract = lambda x: int(x.split('_')[-1].replace('.pickle', ''))
  seed = sorted(list(map(extract, pkl_files)))[0]
  if sep == 'multimer':
    rename = lambda x: f"result_model_{x.split('model_')[1][0]}_multimer_v{x.split('multimer_v')[1][0]}\
_pred_{int(x.split('seed_')[1].split('.')[0]) + pred_shift - seed}.pkl"
  elif sep == 'ptm':
    rename = lambda x: f"result_model_{x.split('model_')[1][0]}_ptm\
_pred_{int(x.split('seed_')[1].split('.')[0]) + pred_shift - seed}.pkl"
  
  new_names = { old: rename(old) for old in pkl_files }
  if FLAGS.do_rename:
    for old in pkl_files:
      new = new_names[old]
      mv(f"{output_path}/{old}", f"{output_path}/{new}")

  return new_names

def rename_pdb(pdb_files:list, output_path:str, pred_shift:int, sep:str):
  extract = lambda x: int(x.split('_')[-1].replace('.pdb', ''))
  seed = sorted(list(map(extract, pdb_files)))[0]
  print(f"Seed used: {seed}")
  if sep == 'multimer':
    rename = lambda x: f"unrelaxed_model_{x.split('model_')[1][0]}\
_multimer_v{x.split('multimer_v')[1][0]}_pred_{int(x.split('seed_')[1].split('.')[0]) + pred_shift - seed}.pdb"
  elif sep == 'ptm':
    rename = lambda x: f"unrelaxed_model_{x.split('model_')[1][0]}\
_ptm_pred_{int(x.split('seed_')[1].split('.')[0]) + pred_shift - seed}.pdb"
    
  new_names = { old: rename(old) for old in pdb_files }
  if FLAGS.do_rename:
    for old in pdb_files:
      new = new_names[old]
      mv(f"{output_path}/{old}", f"{output_path}/{new}")

  map_old_to_new = {
    f"alphafold2{old.split('alphafold2')[1].replace('.pdb', '')}": f"model{new_names[old].split('model')[1].replace('.pdb', '')}"
    for old in new_names
  }

  if FLAGS.do_rename:
    map_file = os.path.join(FLAGS.to_convert, 'unified_map.json')
    if os.path.isfile(map_file):
      with open(map_file, 'r') as map_json:
        name_map = json.load(map_json)
      name_map.update(map_old_to_new)
      with open(map_file, 'w') as map_json:
        json.dump(name_map, map_json, indent=4)
    else:
      with open(map_file, 'w') as map_json:
        json.dump(map_old_to_new, map_json, indent=4)
  return new_names

def create_ranking(predictions_to_rank:pd.core.frame.DataFrame, output_path:str, preset:str):
  for score in ['plddts', 'ptm', 'iptm', 'iptm+ptm']:
    if preset == 'multimer' and score == 'plddts':
      continue
    try:
      df = predictions_to_rank.sort_values(score, ascending=False)
    except KeyError:
      if score == 'plddts':
        print(predictions_to_rank)
        print('No "ptm" found')
        sys.exit()
      continue
    metric_name = score if score != 'iptm+ptm' else 'debug'
    if preset == 'ptm':
      metric_name = 'debug'
    ranking_file_name = f"{output_path}/ranking_{metric_name}.json"
    scores_dict = df.set_index('prediction').to_dict(orient='dict')[score]
    if score != 'iptm+ptm':
      scores_dict = {key: value.item() for key, value in scores_dict.items()}
    order = df['prediction'].to_list()
    
    with open(ranking_file_name, 'w') as json_scores:
      json.dump({ score: scores_dict, 'order': order }, json_scores, indent=4)

def rank_predictions(output_path:str, pdb_files:list, new_pdb_names:list, preset):
  jobname = [ name for name in os.listdir(output_path) if name.endswith('a3m') ]
  jobname = jobname[0].split('.')[0]

  pickle_files = map(lambda x: f"{jobname}_all_{x.replace('.pdb', '.pickle').split('_unrelaxed_')[1]}", pdb_files)
  all_preds = pd.DataFrame()

  for pickle_name, pdb_name in zip(list(pickle_files), new_pdb_names):
    if not os.path.exists(f"{output_path}/{pickle_name}"):
      raise FileNotFoundError(f"{output_path}/{pickle_name}")
    pred_name = f"model_{pdb_name.split('_model_')[1].split('.')[0]}"
    with open(f"{output_path}/{pickle_name}", 'rb') as pickle_scores:
      scores = pickle.load(pickle_scores)
      if preset == 'multimer':
        new_pred = pd.DataFrame([{
          'prediction': pred_name,
          'iptm': scores['iptm'],
          'ptm': scores['ptm'],
          'iptm+ptm': 0.8 * scores['iptm'] + 0.2 * scores['ptm']
          }])
      elif preset == 'ptm':
        new_pred = pd.DataFrame([{
          'prediction': pred_name,
          'plddts': scores['mean_plddt'],
          }])
        
      all_preds = pd.concat([all_preds, new_pred], ignore_index=True)
  
  create_ranking(all_preds, output_path, preset)

def convert_output(output_path:str, pred_shift:int):
  pkls = [ file for file in os.listdir(output_path) if file.endswith('.pickle') ]
  pdbs = [ file for file in os.listdir(output_path) if file.endswith('.pdb') and 'rank' in file ]
  if 'multimer' in pdbs[0]:
    sep = 'multimer'
  elif 'ptm' in pdbs[0]:
    sep = 'ptm'
  else:
    raise ValueError('Neither multimer nor monomer_ptm, an error occured somewhere')

  # rename files
  renamed_pdbs = rename_pdb(pdbs, output_path, pred_shift, sep=sep)
  rank_predictions(output_path, pdbs, renamed_pdbs.values(), preset=sep)
  rename_pkl(pkls, output_path, pred_shift, sep=sep)
    

def move_output(output_path:str, batch):
  whole_path = os.path.realpath(output_path)
  sequence = os.path.basename(os.path.dirname(whole_path))
  batch_path = f"{whole_path}/{batch}"
  destination_path = f"{whole_path}/{batch}/{sequence}"
  to_move = os.listdir(batch_path)
  if FLAGS.do_rename:
    os.makedirs(destination_path)

  for element in to_move:
    source = os.path.join(batch_path, element)
    destination = os.path.join(destination_path, element)
    if FLAGS.do_rename:
      shutil.move(source, destination)

def main(argv):
  if not FLAGS.conversion or not FLAGS.to_convert:
    raise ValueError('Parameter --conversion and --to_convert are mandatory.')
  
  if FLAGS.conversion == 'input':
    assert FLAGS.to_convert.endswith('.fasta') and os.path.isfile(FLAGS.to_convert)
   
    colabfold_input_dir = f"{os.path.dirname(FLAGS.to_convert)}/converted_for_colabfold/" 
    colabfold_fasta_exists = os.path.exists(f"{os.path.dirname(FLAGS.to_convert)}/converted_for_colabfold/")
    if not colabfold_fasta_exists:
      os.makedirs(colabfold_input_dir)

    convert_fasta(FLAGS.to_convert)

  if FLAGS.conversion == 'output':
    assert FLAGS.batches_file
    with open(FLAGS.batches_file, 'r') as batches_json:
      all_batches_infos = json.load(batches_json)

    batches = [ batch for batch in os.listdir(FLAGS.to_convert) if batch.startswith('batch_') ]
    batches = sorted(batches, key=lambda x: int(x.split('_')[1]))

    for batch in batches:
      batch_number = batch.split('_')[1]
      batch_shift = int(all_batches_infos[batch_number]['start'])
      convert_output(f"{FLAGS.to_convert}/{batch}", batch_shift)
      move_output(FLAGS.to_convert, batch)

if __name__ == "__main__": 
  app.run(main)
