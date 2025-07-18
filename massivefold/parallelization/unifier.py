#!/usr/bin/env python
import os
import re
import copy
from absl import flags, app
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
  ['input', 'input_inference', 'output', 'output_singular'],
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
  None,
  ["AFmassive", "ColabFold", "AlphaFold3"],
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
  
  output_fasta = f"{fasta_dir}/converted_for_colabfold/{fasta_file}.fasta"
  with open(output_fasta, 'w') as output:
    output.write(converted)
  
  print(f"Fasta file has been successfully converted for ColabFold at {output_fasta}\n")

def create_alphafold3_json(fasta_path: str, adapted_input_dir: str):
  json_params = os.path.realpath(FLAGS.json_params)
  assert os.path.exists(json_params) and json_params.endswith('.json'), \
  "Please provide a valid path to a json file with --json_params"
  
  all_params = json.load(open(json_params, 'r')) 
  template_dir = all_params['massivefold']['jobfile_templates_dir']
  json_template = os.path.realpath(os.path.join(template_dir, "AlphaFold3", "af3_input.json"))
  json_input = json.load(open(json_template, 'r'))

  parsed_records = list(SeqIO.parse(fasta_path, "fasta"))
  entities = all_params['AF3_run']['fasta_chains']

  # possible entities: "protein", "dna", "rna"
  assert len(parsed_records) == len(entities), \
  f"The number of 'fasta_chains' entities in {json_params} should be the same as in {fasta_path}."
  all_chain_ids = string.ascii_uppercase
  assert len(parsed_records) < len(all_chain_ids), \
  f"Using more than {len(all_chain_ids)} is currently unsupported"

  records = []
  for entity, record in zip(entities, parsed_records):
    records.append({"entity": entity, "sequence_type": "sequence", "seq": str(record.seq)})
  sequence_dicts, bonds = af3_records_to_sequences(records, json_input)

  # create AlphaFold3 input as json
  fasta_file = os.path.basename(fasta_path).split('.fa')[0]
  json_input['name'] = "msas_alphafold3"
  json_input['sequences'] = sequence_dicts
  with open(os.path.join(adapted_input_dir, fasta_file + '.json'), "w") as f:
    json.dump(json_input, f, indent=4)
  print(json.dumps(json_input, indent=4))

def convert_input(args, tool):
  fasta_file = args['input']
  input_dir = os.path.dirname(fasta_file)
  if tool == 'ColabFold':
    adapted_input_dir =  f"{input_dir}/converted_for_colabfold/" 
    os.makedirs(adapted_input_dir, exist_ok=True)
    convert_colabfold_fasta(fasta_file)
  elif tool == 'AlphaFold3':
    adapted_input_dir = f"{input_dir}/alphafold3_json_requests/" 
    os.makedirs(adapted_input_dir, exist_ok=True)
    create_alphafold3_json(fasta_file, adapted_input_dir)

def set_alphafold3_parameters(af3_input: dict, parameters: list):
  for i, sequence in enumerate(af3_input["sequences"]):
    entity = list(sequence.keys())[0]
    if sequence != "protein":
      continue
    for param in parameters:
      if param == "unpairedMsa" or param == "pairedMsa":
        af3_input['sequences'][i][entity][param] = ""
      elif param == "templates":
        af3_input["sequences"][i][entity][param] = []
  return af3_input

def af3_alter_input(batch_input_json, af3_params):
  """Alter the json input (msas and templates) according to the parameters used."""
  alteration_params = ["unpairedMsa", "pairedMsa", "templates"]
  altered_values = {
    "unpairedMsa": "false",
    "pairedMsa": "false",
    "templates": "false"
  }
  alteration = [ i for i in af3_params if i in alteration_params and af3_params[i] == altered_values[i] ]
  if alteration:
    print(f"Input alteration parameters in use: {', '.join(alteration)}")
    batch_input_json = set_alphafold3_parameters(batch_input_json, alteration)
  return batch_input_json


def glycan_traversal(sugar, parent_index, linkage, state):
  """
  Recursively traverse the glycan tree, recording each sugar instance and linkage-based bonds.

  Parameters:
    sugar: glypy.structure.monosaccharide.Monosaccharide
      The current sugar node.
    parent_index: int or None
      The residue index of the parent sugar, if any.
    linkage: glypy.structure.linkage.Linkage or None
      The linkage object connecting the parent to this sugar.
    state: {'ccdCodes': list, 'bondedAtomPairs': list, 'entity': str, 'residue_counter': list}
  """
  # new residue index for the current sugar instance.
  current_index = state['residue_counter'][0]
  state['residue_counter'][0] += 1
  state['ccdCodes'].append(state['map_code'][sugar.serialize(name="iupac_lite")])
  # bond caracterization
  if parent_index is not None and linkage is not None:
    parent_atom, child_atom = "C" + str(linkage.parent_position), "C" + str(linkage.child_position)
    state['bondedAtomPairs'].append([[state['entity'], parent_index, parent_atom], [state['entity'], current_index, child_atom]])
  # continue traversal focusing on childs of current sugar
  for pos, link in sugar.links.items():
    if link.parent.id == sugar.id:
        glycan_traversal(link.child, current_index, link, state)

def af3_resolve_glycan(glycan_str, chain_id):
  from glypy.io import iupac
  parsed_glycan = iupac.loads(glycan_str, dialect="simple")
  iupac_to_ccd = {
    "Gal": "GAL", "a-Gal": "GAL", "b-Gal": "GLB",
    "Glc": "GLC", "a-Glc": "GLC", "b-Glc": "BGC",
    "Man": "MAN", "a-Man": "BMA", "b-Man": "BMA",
    "Fuc": "FUC", "a-Fuc": "FCA", "b-Fuc": "FCB",
    "GlcNAc": "NAG", "Glc2NAc": "NAG",
    "GalNAc": "NGA", "Gal2NAc": "NGA"
  }
  state = {
    'ccdCodes': [],
    'bondedAtomPairs': [],
    'entity': chain_id,
    'residue_counter': [1],
    'map_code': iupac_to_ccd
  }
  if hasattr(parsed_glycan, "root"):
    glycan_traversal(parsed_glycan.root, parent_index=None, linkage=None, state=state)
  else:
    state['ccdCodes'].append(state['map_code'][parsed_glycan.serialize(name="iupac_lite")])

  return state['ccdCodes'], state['bondedAtomPairs']

def af3_records_to_sequences(records, fasta_ids_sequences): 
  glycosylation_attachment = {
    'N': {"atom": 'ND2', "sugar": ["NAG"]},
    'S': {"atom": 'OG', "sugar": ["NGA", "NAG", "MAN"]},
    'T': {"atom": 'OG1', "sugar": ["NGA", "NAG", "MAN"]},
    'K': {"atom": 'O', "sugar": ["GAL"]}
  }

  used_ids = list({}.keys())
  used_ids = list(fasta_ids_sequences.keys())
  sequence_dicts = []
  all_chain_ids = string.ascii_uppercase
  all_chain_ids = [ i for i in all_chain_ids if i not in used_ids ]
  remaining_records = records.copy()
  all_bonds = []
  # Record format: {"entity": entity, "sequence_type": "sequence|ccdCodes|smiles" "seq": record}
  all_sequences = []
  all_entities = []
  for record, chain_id in zip(records, all_chain_ids):
    record_type = record["sequence_type"][:]
    if record["sequence_type"] == "glycosylation":
      chain, position  = used_ids[record["on_chain_index"]], record["at_position"]
      glycosylated_residue = fasta_ids_sequences[chain][position-1]

      ccdCodes, glycan_bondedAtomPairs = af3_resolve_glycan(record["seq"], chain_id) 
      assert glycosylated_residue in glycosylation_attachment, \
      f"'{glycosylated_residue}' is invalid residue for glycosylation. Valid one are among {list(glycosylation_attachment.keys())}"
      assert ccdCodes[0] in glycosylation_attachment[glycosylated_residue]["sugar"], \
      f"Wrong type of glycosylation on {glycosylated_residue}"

      residue_atom = glycosylation_attachment[glycosylated_residue]["atom"]
      if glycan_bondedAtomPairs:
        glycan_root = glycan_bondedAtomPairs[0][0].copy()
        glycan_root[2] = "C1"
      elif not glycan_bondedAtomPairs and len(ccdCodes) == 1:
        glycan_root = [chain_id, 1, "C1"]
      else:
        raise ValueError(f"No bondedAtomPairs found while more than 1 ccdCodes foud: ccdCodes: {ccdCodes}")

      bondedAtomPairs = [[[chain, record["at_position"], residue_atom], glycan_root]]
      bondedAtomPairs.extend(glycan_bondedAtomPairs)
      all_bonds.extend(bondedAtomPairs)
      record_type  = "ccdCodes"
      record["seq"] = ccdCodes

    # either create new entity or add id to existing one
    if record["seq"] not in all_sequences:
      sequence_dicts.append({record["entity"]:  {"id": chain_id, record_type: record["seq"]}})
      all_sequences.append(record["seq"])
      all_entities.append(record["entity"])
    else:
      index = all_sequences.index(record["seq"])
      if isinstance(sequence_dicts[index][all_entities[index]]["id"], str):
        sequence_dicts[index][all_entities[index]]["id"] = [ sequence_dicts[index][all_entities[index]]["id"] ]
      sequence_dicts[index][all_entities[index]]["id"].append(chain_id)
  return sequence_dicts, all_bonds

def af3_entities_to_records(af3_params, fasta_ids_sequences):
  additional_records = []
  ligand = af3_params['ligand']
  all_ptm_types = ["glycosylation", "phosphorylation", "hydroxylation", "methylation", "acetylation","cyclization"]
  modifs_dummy_seq = {"phosphorylation": "PO3", "hydroxylation": "OH", "methylation": "CH3", "acetylation": "CH3CO", "cyclization":"PCA"}

  used_ids = list(fasta_ids_sequences.keys())
  
  for lig in ligand:
    is_ccdcodes = True if lig['ccdCodes'] and lig['ccdCodes'][0] else False
    if is_ccdcodes and lig['smiles']:
      raise ValueError(f"Chose either 'ccdCodes' or 'smiles' for ligand: {lig}")
    elif not is_ccdcodes and not lig["smiles"]:
      continue
    elif is_ccdcodes:
      additional_records.append({"entity": "ligand", "sequence_type": "ccdCodes", "seq": lig["ccdCodes"]})
    elif lig['smiles']:
      additional_records.append({"entity": "ligand", "sequence_type": "smiles", "seq": lig["smiles"]})

  PTMs = af3_params["modifications"]
  if PTMs:
    # parse each chain's list of PTMs
    for i, single_chain_ptms in enumerate(PTMs):
      for ptm in single_chain_ptms:
        record = {}
        record.update({"sequence_type": ptm["type"], "on_chain_index": i})
        # add informations specific to modification
        if ptm["type"] == "glycosylation":
          if not "sequence" in ptm or not ptm["sequence"]:
            print(ptm)
            continue
          record.update({"entity": "ligand", "seq": ptm["sequence"]})

        elif ptm["type"] in all_ptm_types:
          if not ptm["positions"]:
            continue
          ptm["sequence"] = modifs_dummy_seq[ptm["type"]]
          record.update({"entity": "modifications", "seq": ptm["sequence"]})

        # make one record per modification site (when multiple positions specified)
        for pos in list(sorted(map(int, ptm["positions"]))):
          single_record = copy.deepcopy(record)
          single_record["at_position"] = pos
          assert pos <= len(fasta_ids_sequences[used_ids[i]]), \
          f"{ptm['type']} position {pos} is higher than the total length of the sequence ({len(fasta_ids_sequences[used_ids[i]])})"
          additional_records.append(single_record)
        # display modification info
        print(
          f"- {ptm['type'].capitalize()}: \n{ptm['sequence']}"
          f"\nAt positions {', '.join(list(map(str, ptm['positions'])))} on:\n"
          f"{fasta_ids_sequences[used_ids[i]]}"
        )

  PTMs = [ ptm for ptm in additional_records if ptm["sequence_type"] in all_ptm_types ]
  if PTMs:
    print(f"\n{len(PTMs)} PTMs detected")
  else:
    print("No post-translational modification on the sequences.")
  
  return additional_records

def af3_add_modifications(all_modifications, all_sequences):
  """ Add the modifications specified in MassiveFold parameters to the AlphaFold3 input sequences."""
  proteins_modification_codes = {
    'S': {"phosphorylation": 'SEP'},
    'T': {"phosphorylation": 'TPO'},
    'P': {"hydroxylation": "HYP"},
    'R': {"methylation": "MMO"},
    'K': {"acetylation": "ALY", "methylation": "MLZ"},
    'E': {"cyclization": "PCA"},
  }

  dna_modification_codes = {
    "C": {"methylation": "17E"}
  }

  entity_to_modif = {
    "protein": proteins_modification_codes,
    "dna": dna_modification_codes
  }
  # flatten the sequences to have one chain id per entity
  flattened = []
  for sequence in all_sequences:
    entity = list(sequence.keys())[0]
    entity_ids = sequence[entity]["id"]
    if isinstance(entity_ids, list):
      for id in entity_ids:
        entity_elements = {"id": id}
        entity_elements.update({ e: sequence[entity][e] for e in sequence[entity] if e != "id" })
        flattened.append({entity: entity_elements})
    else:
      flattened.append(sequence)
  flattened = sorted(flattened, key=lambda x: list(x.values())[0]['id'])

  # organize modifications to pair one chain with one group of modifs
  modif_per_chain = [ [] for i in range(len(flattened)) ]
  for modif in all_modifications:
    modif_per_chain[modif["on_chain_index"]].append(modif)
  # add modifs to each chains
  modified_sequences = []
  for i, (chain, chain_modifs) in enumerate(zip(flattened, modif_per_chain)):
    entity, entity_modifications = list(chain.keys())[0], []
    chain_data = copy.deepcopy(chain)[entity]
    if entity == "protein":
      key_modif_type, key_modif_position = "ptmType", "ptmPosition"
    elif entity == "rna" or entity == "dna":
      key_modif_type, key_modif_position = "modificationType", "basePosition"
    elif entity == "ligand":
      continue

    for modif in chain_modifs: # add all chain's modifications
      modif_type, modif_position = modif["sequence_type"], modif["at_position"]

      residue_to_modify = chain.copy()[entity]["sequence"][modif_position - 1]
      assert residue_to_modify in entity_to_modif[entity], \
      f"{residue_to_modify} has no supported modifications ({residue_to_modify}{modif_position})"
      assert modif_type in entity_to_modif[entity][residue_to_modify], \
      f"{modif_type} is not a supported modifications for {residue_to_modify}"

      modif_record = {key_modif_type: entity_to_modif[entity][residue_to_modify][modif_type], key_modif_position: modif_position}

      entity_modifications.append(modif_record)
    chain_data["modifications"].extend(entity_modifications)
    modified_sequences.append({entity: chain_data.copy()})
  # flatten sequences: multiple ids becomes one entry per id 
  merged = []
  for i, sequence in enumerate(modified_sequences):
    entity, insert_to = list(sequence.keys())[0], None
    to_add = [sequence[entity]["sequence"], sequence[entity]["modifications"]]
    for n, existing in enumerate(merged):
      compa_entity = list(existing.keys())[0]
      compa_to = [existing[compa_entity]["sequence"], existing[compa_entity]["modifications"]]
      if to_add == compa_to:
        insert_to = n
        break
    if insert_to != None: # same sequence/modifs found
      if isinstance(merged[n][list(merged[n].keys())[0]]["id"], str):
        merged[n][list(merged[n].keys())[0]]["id"] = [merged[n][list(merged[n].keys())[0]]["id"]]
      merged[n][list(merged[n].keys())[0]]["id"].append(sequence[entity]["id"])
    else:
      merged.append(sequence)
  return merged

def af3_sequences_to_ids(batch_input_json):
  fasta_ids_sequences = {}
  map_id_entity = {}
  for chain in batch_input_json["sequences"]:
    if not chain:
      continue
    entity_name = next(iter(chain.keys()))
    entity = chain[entity_name]

    if entity_name in ["protein", "dna", "rna"]:
      seq_key = "sequence"
    elif entity_name == "ligand":
      if "ccdCodes" in entity and entity["ccdCodes"]:
        seq_key = "ccdCodes"
      elif "smiles" in entity and entity["smiles"]:
        seq_key = "smiles"

    light_chain = {entity_name: {i: entity[i] for i in entity if i in ["id", seq_key, "modifications"]}}
    ids = entity["id"]
    if isinstance(ids, list):
      for id in ids:
        fasta_ids_sequences[id] = entity[seq_key]
        map_id_entity[id] = light_chain
    elif isinstance(ids, str):
      try:
        fasta_ids_sequences[ids] = entity[seq_key]
      except KeyError:
        print(entity)
        raise KeyError
      map_id_entity[ids] = light_chain
  return fasta_ids_sequences

def af3_add_input_entity(batch_input_json, af3_params):
  fasta_ids_sequences = af3_sequences_to_ids(batch_input_json)
  additional_records = af3_entities_to_records(af3_params, fasta_ids_sequences)
  # add the sequence modifications
  modifications_types = ["phosphorylation", "methylation", "acetylation", "hydroxylation", "cyclization"]
  new_modifications_records = [ record for record in additional_records if record["sequence_type"] in modifications_types ]
  batch_input_json["sequences"] = af3_add_modifications(new_modifications_records, batch_input_json["sequences"])
  # add the extra sequences (e.g ligands, glycosylation)
  bonds = []
  sequence_types = "ligand"
  new_sequences_records = [ record for record in additional_records if record["entity"] == sequence_types ]
  additional_sequences, bonds = af3_records_to_sequences(new_sequences_records, fasta_ids_sequences)
  if bonds:
    batch_input_json["bondedAtomPairs"] = bonds
  batch_input_json["sequences"].extend(additional_sequences)

  # display the recorded entities for the run launched
  simplified = {"sequences": []}
  for i, obj  in enumerate(batch_input_json["sequences"]):
    entity = list(obj.keys())[0]
    values = {}
    for j in batch_input_json["sequences"][i][entity]:
      if j not in [ "pairedMsa", "unpairedMsa", "templates"]: # skip large entries
        values[j] = batch_input_json["sequences"][i][entity][j]
    simplified["sequences"].append({entity: values})

  if bonds:
    simplified["bondedAtomPairs"] = bonds

  print(json.dumps(simplified, indent=4))
  return batch_input_json


def get_alphafold3_batch_input(input_json: str, params_json: str, batches: str):
  sequence = os.path.basename(os.path.dirname(os.path.dirname(input_json)))

  all_params = json.load(open(params_json, 'r'))
  massivefold_params = all_params["massivefold"]
  af3_params = all_params["AF3_run"]
  logs_dir = massivefold_params["logs_dir"]
  output_dir = massivefold_params["output_dir"]

  # modify input file according to MassiveFold parameters
  batch_input_json = json.load(open(input_json, 'r'))
  batch_input_json = af3_alter_input(batch_input_json, af3_params)
  batch_input_json = af3_add_input_entity(batch_input_json, af3_params)

  # find run_name
  batch_filename = os.path.basename(batches)
  # batch file: <sequence>_<run>_batches.json
  run_name = re.sub(fr"^{sequence}_", "", batch_filename).replace("_batches.json", "")

  # distribute input file for each batches of the run
  batches = json.load(open(batches, 'r'))
  batches_keys = []
  for i in batches:
    batches_keys.extend(list(batches[i].keys()))
  batches_keys = set(batches_keys)

  for batch in batches:
    starting_seed = random.randint(0, 1_000_000)
    num_seeds = int(batches[batch]['end']) - int(batches[batch]['start']) + 1
    model_seeds = [ starting_seed  + i for i in range(num_seeds)]
    single_batch = copy.deepcopy(batch_input_json)
    single_batch['name'] = 'batch_' + batch
    single_batch['modelSeeds'] = model_seeds

    screening_item = {}
    if "smiles" in batches[batch] and batches[batch]["smiles"]:
      screening_item = {"entity": "ligand", "sequence_type": "smiles", "seq": batches[batch]["smiles"]}
    elif "ccdcode" in batches[batch] and batches[batch]["ccdcode"]:
      screening_item = {"entity": "ligand", "sequence_type": "ccdCodes", "seq": batches[batch]["ccdcode"]}
    if screening_item:
      fasta_ids_sequences = af3_sequences_to_ids(single_batch)
      additional_sequences, _= af3_records_to_sequences([screening_item], fasta_ids_sequences)
      single_batch["name"] = batches[batch]["id"]
      single_batch["sequences"].extend(additional_sequences)

    alphafold3_input = os.path.join(output_dir, sequence, run_name, f"af3_batch_{batch}.json")
    json.dump(single_batch, open(alphafold3_input, 'w'), indent=4)

def prepare_inference(args, tool):
  input = args['input']
  params = args['params']
  batches = args['batches']
  if tool == "AlphaFold3":
    get_alphafold3_batch_input(input, params, batches)

def rename_colabfold_pkl(pkl_files:list, output_path:str, pred_shift:int, sep:str):
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
      #mv(f"{output_path}/{old}", f"{output_path}/{new}")
      cp(f"{output_path}/{old}", f"{output_path}/{new}")

  return new_names

def rename_colabfold_pdb(pdb_files:list, output_path:str, pred_shift:int, sep:str):
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
      #mv(f"{output_path}/{old}", f"{output_path}/{new}")
      cp(f"{output_path}/{old}", f"{output_path}/{new}")

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

def rank_colabfold_predictions(output_path:str, pdb_files:list, new_pdb_names:list, preset):
  jobname = [ name for name in os.listdir(output_path) if name.endswith('a3m') ]
  jobname = jobname[0].split('.')[0]

  pickle_files = map(lambda x: x.replace("_unrelaxed_", "_all_").replace('.pdb', '.pickle'), pdb_files)
  all_preds = pd.DataFrame()

  for pickle_name, pdb_name in zip(list(pickle_files), new_pdb_names):
    if not os.path.exists(f"{output_path}/{pickle_name}"):
      raise FileNotFoundError(f"{output_path}/{pickle_name}")
    pred_name = f"model_{pdb_name.split('_model_')[1].split('.')[0]}"
    with open(f"{output_path}/{pickle_name}", 'rb') as pickle_scores:
      scores = pickle.load(pickle_scores)
    content = {"prediction": pred_name, "ptm": scores["ptm"]}
    if "iptm" in scores:
      content.update({"iptm": scores["iptm"], 'iptm+ptm': 0.8 * scores['iptm'] + 0.2 * scores['ptm']})
    if "mean_plddt" in scores:
      content.update({"plddts": scores["mean_plddt"]})
    new_pred = pd.DataFrame([content])
    all_preds = pd.concat([all_preds, new_pred], ignore_index=True)
  
  create_colabfold_ranking(all_preds, output_path, preset)

def convert_output(tool):
  if tool not in ["ColabFold", "AlphaFold3"]:
    print(f"No conversion needed for {tool} output")
    return
  with open(FLAGS.batches_file, 'r') as batches_json:
    all_batches_infos = json.load(batches_json)

  if tool == "AlphaFold3":
    batches_files = [
      batch_file for batch_file in os.listdir(FLAGS.to_convert)
      if batch_file.startswith('af3_batch_') and batch_file.endswith('.json')
    ]
    batches_files = sorted(
      batches_files,
      key=lambda x: int(x.replace('af3_batch_', '').replace('.json', ''))
    )
    batches = [
      json.load(open(os.path.join(FLAGS.to_convert, batch_file), 'r'))['name']
      for batch_file in batches_files
    ]
  elif tool in ["AFmassive", "ColabFold"]:
    batches = [ batch for batch in os.listdir(FLAGS.to_convert) if batch.startswith('batch_') ]
    batches = sorted(batches, key=lambda x: int(x.split('_')[1]))
    batches_files = batches

  working, not_working = [], []
  for file, batch in zip(batches_files, batches):
    if tool == "ColabFold":
      batch_number = batch.split('_')[1]
      batch_shift = int(all_batches_infos[batch_number]['start'])
      convert_colabfold_output(f"{FLAGS.to_convert}/{batch}", batch_shift)
      move_output(FLAGS.to_convert, batch)
    elif tool == "AlphaFold3":
      batch_number = file.replace('af3_batch_', '').replace('.json', '')
      batch_shift = int(all_batches_infos[batch_number]['start'])
      try:
        convert_alphafold3_output(f"{FLAGS.to_convert}/{batch}", batch_shift)
        working.append(batch)
      except FileNotFoundError:
        not_working.append(batch)

  if not_working:
    print(f"Batch not completed: {' - '.join(not_working)}")


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
  renamed_pdbs = rename_colabfold_pdb(pdbs, output_path, pred_shift, sep=sep)
  rank_colabfold_predictions(output_path, pdbs, renamed_pdbs.values(), preset=sep)
  rename_colabfold_pkl(pkls, output_path, pred_shift, sep=sep)

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
  num_samples = 1 + max([ int(os.path.basename(x).split('sample-')[1]) for x in seed_dirs ])

  seeds_to_0 = [ seed - seeds[0] for seed in seeds ]
  pred_nb = [ pred*num_samples for pred in seeds_to_0 ]
  samples = [ i%num_samples for i, _ in enumerate(pred_nb) ]
  pred_nb = [ pred + i%num_samples for i, pred in enumerate(pred_nb) ]
  
  df = pd.DataFrame( {
    "pred_nb": pred_nb,
    "original_dir": seed_dirs,
    "seed": seeds,
    "sample": samples 
  })
  df = df.merge(df_ranking_scores, on=["seed", "sample"], how="left")
  df["pred_nb"] += pred_shift * num_samples
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
  df = af3_extract_plddts_create_pkl(df, output_path)
  all_score_types = [ i for i in ['ranking_score', 'iptm', 'ptm', 'mean_plddt'] if i in df.columns ]
  df = df.sort_values(all_score_types, ascending=False, ignore_index=True)
  df['rank'] = df.index
  df["ranked_name"] = "ranked_" + df["rank"].astype(str) + "_" + df["prediction_name"] + ".cif"
  af3_move_and_rename(df, output_path)

def prediction_metrics(input_dir: str, nature: str):
  prediction_confs = os.path.join(input_dir, 'summary_confidences.json')
  with open(prediction_confs, 'r') as f:
    data = json.load(f)
  metrics_possibilities = [ "plddt", "iptm", "ptm" ]#, "chain_pair_iptm" ]
  metrics = { metric: data[metric] for metric in metrics_possibilities if metric in data and data[metric]}
  return metrics

def plddts_from_cif(cif_filename):
  from Bio.PDB.MMCIFParser import MMCIFParser
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
  score_types = [ i for i in score_types if i in df.columns ]
  all_scores = { score_map[stype]: { stype: {}, "order": [] } for stype in score_types if stype in df.columns  }

  for pred in pred_list:

    model_cif_name = os.path.join(pred["original_dir"], 'model.cif')
    new_cif_name = os.path.join(os.path.dirname(pred["original_dir"]), pred["ranked_name"])
    cp(model_cif_name, new_cif_name)

    path_to_confidence = os.path.join(os.path.dirname(pred["original_dir"]), "confidences")
    if not os.path.exists(path_to_confidence):
      os.makedirs(path_to_confidence)
    model_confidence_file = os.path.join(pred["original_dir"], 'summary_confidences.json')
    new_confidence_file = os.path.join(os.path.dirname(pred["original_dir"]), "confidences", f'{pred["prediction_name"]}.json')
    cp(model_confidence_file, new_confidence_file)

    model_no_rank = os.path.join(os.path.dirname(pred["original_dir"]), '_'.join(os.path.basename(new_cif_name).split('_')[2:]))
    cp(new_cif_name, model_no_rank)
    for stype in score_types:
      all_scores[score_map[stype]][stype][pred["prediction_name"]] = pred[stype]
      all_scores[score_map[stype]]["order"].append(pred["prediction_name"])

  for score_ranking in all_scores:
    ranking_file = os.path.join(output_dir, f"ranking_{score_ranking}.json")
    json.dump(all_scores[score_ranking], open(ranking_file, 'w'), indent=4)

def af3_extract_plddts_create_pkl(df, output_dir):
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

def main(argv):
  assert FLAGS.conversion and FLAGS.to_convert, \
  'Parameter --conversion and --to_convert are mandatory.'
  assert FLAGS.tool, "Please specify the tool used for prediction."
  
  if FLAGS.conversion == 'input':
    tools = ["ColabFold", "AlphaFold3"]
    if FLAGS.tool not in tools:
      print(f"No input conversion needed with {FLAGS.tool}.")
      return 
    assert os.path.isfile(FLAGS.to_convert) and FLAGS.to_convert.endswith('.fasta'), \
    f"Fasta file is invalid {FLAGS.to_convert}"
    convert_input({"input": FLAGS.to_convert}, FLAGS.tool)

  elif FLAGS.conversion == 'input_inference':
    tools = ["AlphaFold3"]
    if FLAGS.tool not in tools:
      print(f"No input preparation in anticipation of inference for {FLAGS.tool}.")
      return  
    assert os.path.isfile(FLAGS.batches_file) and FLAGS.batches_file.endswith('.json'), \
    f"Json batches file is invalid: {FLAGS.batches_file}"
    assert os.path.isfile(FLAGS.to_convert) and FLAGS.to_convert.endswith('.json'), \
    f"Json request file is invalid: {FLAGS.to_convert}"
    arguments = {"input": FLAGS.to_convert, "params": FLAGS.json_params, "batches": FLAGS.batches_file}
    prepare_inference(arguments, FLAGS.tool)

  elif FLAGS.conversion == 'output':
    tools = ["ColabFold", "AlphaFold3"]
    if FLAGS.tool not in tools:
      print(f"No output standardization for {FLAGS.tool}.")
      return 
    assert FLAGS.batches_file, 'Json batches file (--batches_file) is mandatory for output conversion (--conversion output)'
    print(f"Convert for tool {FLAGS.tool}")
    convert_output(FLAGS.tool)

  elif FLAGS.conversion == "output_singular":
    tools = ["ColabFold", "AlphaFold3"]
    if FLAGS.tool not in tools:
      print(f"No output standardization for {FLAGS.tool}.")
      return
    if FLAGS.tool == "ColabFold":
      convert_colabfold_output(FLAGS.to_convert, 0)
      move_output(os.path.dirname(os.path.realpath(FLAGS.to_convert)), "batch_0")
    elif FLAGS.tool == "AlphaFold3":
      convert_alphafold3_output(FLAGS.to_convert, 0)

if __name__ == "__main__": 
  app.run(main)
