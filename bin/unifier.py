#!/usr/bin/env python
import os
from absl import flags, app
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import SeqIO
import pickle
from shutil import move as mv, copy as cp
import pandas as pd
import json
import sys
import shutil
import string
import random
import numpy as np

FLAGS = flags.FLAGS
flags.DEFINE_enum(
  'conversion',
  None,
  ['input', 'input_inference', 'output'],
  "What to convert."
  "'input' to get the fasta format of colabfold from a traditionnal multichain pdb."
  "'output' to transform colabfold output to the alphafold one.")
flags.DEFINE_string(
  'to_convert',
  '',
  'Path of the fasta file (for --conversion=input) or output directory (for --conversion=output)'
  'To convert')
flags.DEFINE_enum(
  'tool',
  "ColabFold",
  ["ColabFold", "AlphaFold3"],
  "Chose the tool from which the input/output should be unified.")
flags.DEFINE_string(
  "json_params",
  "",
  "Set json file path for input parameters. "
  "Necessary when using '--tool AlphaFold3' coupled with '--conversion input'")
flags.DEFINE_string(
  'batches_file',
  '',
  'Path to batches file. If --conversion=output, this file is necessary.')
flags.DEFINE_bool(
  'do_rename',
  True,
  'To rename file or not')

def convert_colabfold_fasta(fasta_path:str):
  records = list(SeqIO.parse(fasta_path, "fasta"))
  fasta_dir = os.path.dirname(fasta_path)
  fasta_file = os.path.basename(fasta_path).split('.fa')[0]

  converted = ''
  ids = '>'
  sequences = ''
  for record in records:
    ids += f"{record.id}:"
    sequences += f"{record.seq}:"

  converted = f"{ids[:-1]}\n{sequences[:-1]}"
  
  with open(f"{fasta_dir}/converted_for_colabfold/{fasta_file}.fasta", 'w') as output:
    output.write(converted)

def create_alphafold3_json(fasta_path: str, adapted_input_dir: str):
  json_params = os.path.realpath(FLAGS.json_params)
  assert os.path.exists(json_params) and json_params.endswith('.json'), \
  "Please provide a valid path to a json file with --json_params"
  
  template_dir = json.load(open(json_params, 'r'))['massivefold']['jobfile_templates_dir']
  json_template = os.path.realpath(os.path.join(template_dir, "AlphaFold3", "af3_input.json"))
  with open(json_template, 'r') as template:
    input_json = json.load(template)

  all_chain_ids = string.ascii_uppercase + string.ascii_lowercase
  records = list(SeqIO.parse(fasta_path, "fasta"))
  entities = json.load(open(json_params, 'r'))['AF3_run']['entities']

  assert len(records) == len(entities), \
  f"The number of entities in {json_params} should be the same as in {fasta_path}."
  assert len(records) < len(all_chain_ids), \
  f"Using more than {len(all_chain_ids)} is currently unsupported"
  
  sequence_dicts = []
  for entity, record, chain_id in zip(entities, records, all_chain_ids):
    all_sequences = [sequence[list(sequence.keys())[0]]["sequence"] for sequence in sequence_dicts]
    if f"{record.seq}" not in all_sequences:
      sequence_dicts.append({entity:  {"id": [chain_id], "sequence": f"{record.seq}"}})
    else:
      index = all_sequences.index(f"{record.seq}")
      sequence_dicts[index][entities[index]]["id"].append(chain_id)

  # create AlphaFold3 input as json
  fasta_file = os.path.basename(fasta_path).split('.fa')[0]
  json_input = json.load(open(json_template, 'r'))
  json_input['name'] = "msas_alphafold3"
  json_input['sequences'] = sequence_dicts
  with open(os.path.join(adapted_input_dir, fasta_file + '.json'), "w") as f:
    json.dump(json_input, f, indent=4)

def convert_input(args, tool):
  fasta_file = args['input']
  input_dir = os.path.dirname(fasta_file)

  if tool == 'ColabFold':
    adapted_input_dir =  f"{input_dir}/converted_for_colabfold/" 
    if not os.path.exists(adapted_input_dir):
      os.makedirs(adapted_input_dir)
    convert_colabfold_fasta(fasta_file)

  elif tool == 'AlphaFold3':
    adapted_input_dir = f"{input_dir}/alphafold3_json_requests/" 
    if not os.path.exists(adapted_input_dir):
      os.makedirs(adapted_input_dir)
    create_alphafold3_json(fasta_file, adapted_input_dir)

def get_alphafold3_batch_input(input_json: str, params_json: str, batches: str):
  sequence = os.path.basename(os.path.dirname(os.path.dirname(input_json)))
  logs_dir = json.load(open(params_json, 'r'))["massivefold"]["logs_dir"]
  output_dir = json.load(open(params_json, 'r'))["massivefold"]["output_dir"]
  batch_filename = os.path.basename(batches)
  run_name = batch_filename.replace(f"{sequence}_", "").replace("_batches.json", "")

  batches = json.load(open(batches, 'r'))
  for batch in batches:
    starting_seed = random.randint(0, 1_000_000)
    num_seeds = int(batches[batch]['end']) - int(batches[batch]['start']) + 1
    model_seeds = [ starting_seed  + i for i in range(num_seeds)]
    alphafold3_input = os.path.join(
      output_dir,
      sequence,
      run_name,
      f"af3_batch_{batch}.json")
    #os.makedirs(os.path.dirname(alphafold3_input))
    batch_input_json = json.load(open(input_json, 'r'))
    batch_input_json['name'] = 'batch_' + batch
    batch_input_json['modelSeeds'] = model_seeds
    json.dump(batch_input_json, open(alphafold3_input, 'w'), indent=4)

def prepare_inference(args, tool):
  input = args['input']
  params = args['params']
  batches = args['batches']
  if tool == "AlphaFold3":
    get_alphafold3_batch_input(input, params, batches)


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

def create_colabfold_ranking(predictions_to_rank:pd.core.frame.DataFrame, output_path:str, preset:str):
  metrics = ['ptm', 'iptm', 'iptm+ptm'] if preset == 'multimer' else ['plddts', 'ptm']
  
  for metric in metrics:
    try:
      df = predictions_to_rank.sort_values(metric, ascending=False)
    except KeyError:
      print(predictions_to_rank)
      print(f'No "{metric}" found')
      sys.exit()

    if preset == 'multimer':
      metric_name = metric if metric != 'iptm+ptm' else 'debug'
    elif preset == 'ptm':
      metric_name = metric if metric != 'plddts' else 'debug'

    ranking_file_name = f"{output_path}/ranking_{metric_name}.json"
    scores_dict = df.set_index('prediction').to_dict(orient='dict')[metric]
    if metric != 'iptm+ptm':
      scores_dict = {key: value.item() for key, value in scores_dict.items()}
    order = df['prediction'].to_list()
    
    with open(ranking_file_name, 'w') as json_scores:
      json.dump({ metric: scores_dict, 'order': order }, json_scores, indent=4)

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
          'ptm': scores['ptm'],
          'plddts': scores['mean_plddt'],
          }])
        
      all_preds = pd.concat([all_preds, new_pred], ignore_index=True)
  
  create_colabfold_ranking(all_preds, output_path, preset)

def convert_colabfold_output(output_path:str, pred_shift:int):
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

def prediction_metrics(input_dir: str, nature: str):
  prediction_confs = os.path.join(input_dir, 'summary_confidences.json')
  with open(prediction_confs, 'r') as f:
    data = json.load(f)
  metrics_possibilities = [ "plddt", "iptm", "ptm" ]#, "chain_pair_iptm" ]
  metrics = { metric: data[metric] for metric in metrics_possibilities if metric in data }
  return metrics

def plddts_from_cif(cif_filename):
  pdb_parser = MMCIFParser(QUIET=True)
  structure = pdb_parser.get_structure(
    'strct',
    cif_filename)
  model = next(structure.get_models())
  residues = model.get_residues()
  residues_plddt = []
  for res in residues:
    residues_plddt.append(np.round(np.mean([ i.get_bfactor() for i in res.get_atoms() ]), decimals=2))
  return residues_plddt

def af3_move_and_rename(df, output_dir):
  pred_list = df.to_dict(orient="records")
  score_map = { "ranking_score": "debug", "iptm": "iptm", "ptm": "ptm", "mean_plddt": "plddt"  }
  score_types = ['iptm', 'ptm', 'ranking_score', 'mean_plddt' ]
  all_scores = { score_map[stype]: { stype: {}, "order": [] } for stype in score_types if stype in df.columns  }

  for pred in pred_list:
    model_cif_name = os.path.join(pred["original_dir"], 'model.cif')
    new_cif_name = os.path.join(os.path.dirname(pred["original_dir"]), pred["ranked_name"])
    cp(model_cif_name, new_cif_name)
    model_no_rank = os.path.join(os.path.dirname(pred["original_dir"]), '_'.join(os.path.basename(new_cif_name).split('_')[2:]))
    cp(new_cif_name, model_no_rank)
    for stype in score_types:
      all_scores[score_map[stype]][stype][pred["prediction_name"]] = pred[stype]
      all_scores[score_map[stype]]["order"].append(pred["prediction_name"])

  for score_ranking in all_scores:
    ranking_file = os.path.join(output_dir, f"ranking_{score_ranking}.json")
    json.dump(all_scores[score_ranking], open(ranking_file, 'w'), indent=4)

def exctract_plddts_create_pkl(df, output_dir):
  pred_list = df.to_dict(orient="records")
  for pred in pred_list:
    model_cif_name = os.path.join(pred["original_dir"], 'model.cif')
    pred_plddts = plddts_from_cif(model_cif_name)
    pred["mean_plddt"] = np.mean(pred_plddts)

    json_confidences_name = os.path.join(pred["original_dir"], 'confidences.json')
    json_confidences_file = json.load(open(json_confidences_name, 'r'))
    json_confidences_file["predicted_aligned_error"] = json_confidences_file["pae"]
    json_confidences_file["max_predicted_aligned_error"] = np.max(json_confidences_file["pae"])
    json_confidences_file["plddt"] = pred_plddts

    pkl_name = os.path.join(os.path.dirname(pred["original_dir"]), pred["pkl_name"])
    pickle.dump(json_confidences_file, open(pkl_name, 'wb'))

  updated_df = pd.DataFrame(pred_list)
  return updated_df

def convert_alphafold3_output(output_path: str, pred_shift: int):
  df_ranking_scores = pd.read_csv(os.path.join(output_path, "ranking_scores.csv"))
  seed_dirs = [ os.path.join(output_path, seed) for seed in os.listdir(output_path) if seed.startswith('seed-')]
  seed_dirs = sorted(
    seed_dirs,
    key=lambda x: (
      int(os.path.basename(x).replace('seed-', '').split('_sample')[0]),
      int(os.path.basename(x).split('_sample-')[1])
    )
  )
  seeds = [ int(os.path.basename(x).replace('seed-', '').split('_sample')[0]) for x in seed_dirs ]
  seeds_to_0 = [ seed - seeds[0] for seed in seeds ]
  pred_nb = [ pred*5 for pred in seeds_to_0 ]
  samples = [ i%5 for i, _ in enumerate(pred_nb) ]
  pred_nb = [ pred + i%5 for i, pred in enumerate(pred_nb) ]
  
  df = pd.DataFrame( {
    "pred_nb": pred_nb,
    "original_dir": seed_dirs,
    "seed": seeds,
    "sample": samples 
  })
  df = df.merge(df_ranking_scores, on=["seed", "sample"], how="left")
  df["pred_nb"] += pred_shift * 5
  to_add = {}
  for seed in seed_dirs:
    pred_metrics = prediction_metrics(seed, "multimer")
    for metric in pred_metrics:
      if metric not in to_add:
        to_add[metric] = [pred_metrics[metric]]
      else:
        to_add[metric].append(pred_metrics[metric])
  for new_metric in to_add:
    df[new_metric] = to_add[new_metric]

  df["prediction_name"] = "af3" + "_seed_" + df["seed"].astype(str) \
    + "_sample_" + df["sample"].astype(str) + "_pred_" + df["pred_nb"].astype(str)
  df["pkl_name"] = "result_" + df["prediction_name"] + ".pkl"
  df = exctract_plddts_create_pkl(df, output_path)
  df = df.sort_values(['ranking_score', 'iptm', 'ptm', 'mean_plddt'], ascending=False, ignore_index=True)
  df['rank'] = df.index
  df["ranked_name"] = "ranked_" + df["rank"].astype(str) + "_" + df["prediction_name"] + ".cif"
  af3_move_and_rename(df, output_path)

def main(argv):
  
  assert FLAGS.conversion and FLAGS.to_convert, \
  'Parameter --conversion and --to_convert are mandatory.'
  
  if FLAGS.conversion == 'input':
    assert os.path.isfile(FLAGS.to_convert) and FLAGS.to_convert.endswith('.fasta')
    convert_input({"input": FLAGS.to_convert}, FLAGS.tool)

  elif FLAGS.conversion == 'input_inference':
    assert os.path.isfile(FLAGS.batches_file) and FLAGS.batches_file.endswith('.json')
    assert os.path.isfile(FLAGS.to_convert) and FLAGS.to_convert.endswith('.json')
    prepare_inference({"input": FLAGS.to_convert,
                       "params": FLAGS.json_params,
                       "batches": FLAGS.batches_file},
                      FLAGS.tool)

  elif FLAGS.conversion == 'output':
    assert FLAGS.batches_file
    with open(FLAGS.batches_file, 'r') as batches_json:
      all_batches_infos = json.load(batches_json)

    batches = [ batch for batch in os.listdir(FLAGS.to_convert) if batch.startswith('batch_') ]
    batches = sorted(batches, key=lambda x: int(x.split('_')[1]))

    for batch in batches:
      batch_number = batch.split('_')[1]
      batch_shift = int(all_batches_infos[batch_number]['start'])
      if FLAGS.tool == "ColabFold":
        convert_colabfold_output(f"{FLAGS.to_convert}/{batch}", batch_shift)
        move_output(FLAGS.to_convert, batch)
      elif FLAGS.tool == "AlphaFold3":
        convert_alphafold3_output(f"{FLAGS.to_convert}/{batch}", batch_shift)

if __name__ == "__main__": 
  app.run(main)
