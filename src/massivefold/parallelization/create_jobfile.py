#!/usr/bin/env python

import argparse
import json
from string import Template
import shutil
import os

def create_single_jobfile(
  jobfile_type,
  templates,
  params,
  json_params,
  mf_following_msas,
  mf_before_inference,
  create_files,
  sequence_name,
  run_name):
  params["json_params"] = json_params

  if jobfile_type == "alignment":
    params["mf_following_msas"] = mf_following_msas
  elif jobfile_type == "jobarray":
    params["mf_before_inference"] = mf_before_inference

  jobfile = Template(templates[jobfile_type]).substitute(params)
  if create_files:
    with open(f"{sequence_name}_{run_name}_{jobfile_type}.slurm", 'w') as slurm_job:
      slurm_job.write(jobfile)
  return jobfile

def create_all_jobfile(
  templates,
  params,
  json_params,
  mf_following_msas,
  mf_before_inference,
  create_files,
  sequence_name,
  run_name):
  for jobtype in ['alignment', 'jobarray', 'post_treatment']:
    create_single_jobfile(
      jobtype,
      templates,
      params,
      json_params,
      mf_following_msas,
      mf_before_inference,
      create_files,
      sequence_name,
      run_name)

def group_templates(all_params, job_types, tool):
  templates_paths = all_params['massivefold']
  sequence = all_params['massivefold']['sequence_name']
  run = all_params['massivefold']['run_name']
  grouped_templates = {}

  if tool == "AFmassive":
    tool_code = "AFM"
    model_preset = all_params[f'{tool_code}_run'][f'model_preset']
  elif tool == "AlphaFold3":
    tool_code = "AF3"
    model_preset = "multimer"
  elif tool == "ColabFold":
    tool_code = "CF"
    model_preset = all_params[f'{tool_code}_run'][f'model_preset']

  for job_type in job_types:
    template_file = f"{job_type}_{model_preset}"
    header_path = f"{templates_paths['jobfile_headers_dir']}/{job_type}.slurm"
    template_path = f"{templates_paths['jobfile_templates_dir']}/{tool}/{template_file}.slurm"
    jobfile = f"{sequence}_{run}_{job_type}.slurm"
    merge_header_and_template(header_path, template_path, jobfile)

    with open(jobfile, 'r') as temp_file:
      jobfile = temp_file.read() 
    grouped_templates[job_type] = jobfile
  return grouped_templates

def merge_header_and_template(header, template, jobfile_name):
  with open(jobfile_name, 'wb') as outfile:
    for filename in [header, template]:
      with open(filename, 'rb') as infile:
        shutil.copyfileobj(infile, outfile)

def detect_model_preset(fasta_file):
  with open(fasta_file, 'r') as fasta_content:
    lines = fasta_content.readlines()

  model_preset = "monomer_ptm"
  sequence_number = 0
  for line in lines:
    if line.startswith('>'):
      sequence_number += 1
    if sequence_number > 1:
      model_preset = "multimer"
      break
  return model_preset

def main(
  job_type,
  sequence_name,
  run_name,
  mf_following_msas,
  mf_before_inference,
  create_files,
  path_to_parameters,
  tool):
  # extract the number of batch for the jobarray
  with open(f'{sequence_name}_{run_name}_batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
  run_params = {
      'run_name': run_name,
      'sequence_name': sequence_name,
      'substitute_batch_number': list(batches.keys())[-1]
    }

  # parameters parsing, computing and display
  if not path_to_parameters:
    raise ValueError('Parameters files missing, use --path_to_parameters.')

  if tool == "AFmassive":
    tool_code = "AFM"
  elif tool == "AlphaFold3":
    tool_code = "AF3"
  elif tool == "ColabFold":
    tool_code = "CF"
  
  with open(path_to_parameters, 'r') as parameters_json:
    all_params = json.load(parameters_json)

  all_params['massivefold']['sequence_name'] = sequence_name
  all_params['massivefold']['run_name'] = run_name

  preset_dict = {"model_preset": detect_model_preset(
    os.path.join(all_params["massivefold"]["input_dir"], f"{sequence_name}.fasta")
  )}

  all_params[f"{tool_code}_run"] = preset_dict | all_params[f"{tool_code}_run"]
  run_params.update(all_params['massivefold'])
  run_params.update(all_params['custom_params'])
  run_params.update(all_params[f'{tool_code}_run'])
  run_params.update(all_params['plots'])

  if job_type == "jobarray":
    print("Parameters of the run:")
    for i in all_params['custom_params']:
      print(f"{i}: {all_params['custom_params'][i]}")
    for i in all_params[f'{tool_code}_run']:
      print(f"{i}: {all_params[f'{tool_code}_run'][i]}")
    print()
  if job_type == "post_treatment":
    print("Parameters for plots:")
    for i in all_params['plots']:
      print(f"{i}: {all_params['plots'][i]}")

  if job_type != 'all':
    all_templates = group_templates(all_params, [job_type], tool)
    create_single_jobfile(
      job_type,
      all_templates,
      run_params,
      path_to_parameters,
      mf_following_msas,
      mf_before_inference,
      create_files,
      sequence_name,
      run_name)
  else:
    all_templates = group_templates(all_params, ['alignment', 'jobarray', 'post_treatment'], tool)
    create_all_jobfile(
      all_templates,
      run_params,
      path_to_parameters,
      mf_following_msas,
      mf_before_inference,
      create_files,
      sequence_name,
      run_name)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--job_type',
    default='all',
    choices=['all', 'alignment', 'jobarray', 'post_treatment'],
    help='Type of the jobfile to create')
  parser.add_argument('--sequence_name', default='', help='name of the fasta sequence used for the run.')
  parser.add_argument(
    '--run_name',
    default='',
    help='name of the run, it can be anything and it will be the name of the output directory under the sequence name.')
  parser.add_argument(
    '--mf_following_msas',
    default='true',
    choices=['true', 'false'],
    help='to activate if not using only_msas params.')
  parser.add_argument(
    '--mf_before_inference',
    default='false',
    choices=['true', 'false'],
    help='to activate if using jobid (-j) params.')
  parser.add_argument(
    '--create_files',
    action=argparse.BooleanOptionalAction,
    default=True,
    help='')
  parser.add_argument(
    '--path_to_parameters',
    default='',
    help='Path to a json file were the jobfile parameters can be specified.')
  parser.add_argument(
    '--tool',
    required=True,
    choices=['AFmassive', 'ColabFold', 'AlphaFold3'],
    help='Specify the tool used by MassiveFold for structure prediction.')

  parsed = parser.parse_args()
  main(
    parsed.job_type,
    parsed.sequence_name,
    parsed.run_name,
    parsed.mf_following_msas,
    parsed.mf_before_inference,
    parsed.create_files,
    parsed.path_to_parameters,
    parsed.tool
  )
