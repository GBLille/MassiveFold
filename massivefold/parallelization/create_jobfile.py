#!/usr/bin/env python

import json
from string import Template
from absl import app, flags
import shutil

FLAGS = flags.FLAGS
flags.DEFINE_enum(
  "job_type",
  'all',
  ['all', 'alignment', 'jobarray', 'post_treatment'],
  'Type of the jobfile to create')
flags.DEFINE_string(
  "sequence_name",
  "",
  'name of the fasta sequence used for the run.')
flags.DEFINE_string(
  "run_name",
  '',
  'name of the run, it can be anything and it will'
  'be the name of the output directory under the sequence name.')
flags.DEFINE_bool("create_files", True, '')
flags.DEFINE_string(
  "path_to_parameters",
  "",
  "Path to a json file were the jobfile parameters can be specified.")
flags.DEFINE_enum(
  "tool",
  "AFmassive",
  ["AFmassive", "ColabFold", "AlphaFold3"],
  "Specify the tool used by MassiveFold for structure prediction.")

def create_single_jobfile(jobfile_type, templates:dict, params, json_params: str):
  params["json_params"] = json_params
  jobfile = Template(templates[jobfile_type]).substitute(params)
  if FLAGS.create_files:
    with open(f"{FLAGS.sequence_name}_{FLAGS.run_name}_{jobfile_type}.slurm", 'w') as slurm_job:
      slurm_job.write(jobfile)
  return jobfile

def create_all_jobfile(templates:dict, params:dict, json_params: str):
  for jobtype in ['alignment', 'jobarray', 'post_treatment']:
    create_single_jobfile(jobtype, templates, params, json_params)

def group_templates(all_params, job_types:list):
  templates_paths = all_params['massivefold']
  sequence = all_params['massivefold']['sequence_name']
  run = all_params['massivefold']['run_name']
  grouped_templates = {}
  
  if FLAGS.tool == "AFmassive":
    tool_code = "AFM"
  elif FLAGS.tool == "AlphaFold3":
    tool_code = "AF3"
  elif FLAGS.tool == "ColabFold":
    tool_code = "CF"
  
  model_preset = all_params[f'{tool_code}_run'][f'model_preset']
  for job_type in job_types:
    template_file = f"{job_type}_{model_preset}"
    header_path = f"{templates_paths['jobfile_headers_dir']}/{job_type}.slurm"
    template_path = f"{templates_paths['jobfile_templates_dir']}/{FLAGS.tool}/{template_file}.slurm"
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

def main(argv):
  # extract the number of batch for the jobarray
  with open(f'{FLAGS.sequence_name}_{FLAGS.run_name}_batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
  run_params = {
      'run_name': FLAGS.run_name,
      'sequence_name': FLAGS.sequence_name,
      'substitute_batch_number': list(batches.keys())[-1]
    }

  # parameters parsing, computing and display
  if not FLAGS.path_to_parameters:
    raise ValueError('Parameters files missing, use --path_to_parameters.')

  if FLAGS.tool == "AFmassive":
    tool_code = "AFM"
  elif FLAGS.tool == "AlphaFold3":
    tool_code = "AF3"
  elif FLAGS.tool == "ColabFold":
    tool_code = "CF"
  
  with open(FLAGS.path_to_parameters, 'r') as parameters_json:
    all_params = json.load(parameters_json)

  all_params['massivefold']['sequence_name'] = FLAGS.sequence_name
  all_params['massivefold']['run_name'] = FLAGS.run_name

  run_params.update(all_params['massivefold'])
  run_params.update(all_params['custom_params'])
  run_params.update(all_params[f'{tool_code}_run'])
  run_params.update(all_params['plots'])

  if FLAGS.job_type == "jobarray":  
    print("Parameters of the run:")
    for i in all_params['custom_params']:
      print(f"{i}: {all_params['custom_params'][i]}")
    for i in all_params[f'{tool_code}_run']:
      print(f"{i}: {all_params[f'{tool_code}_run'][i]}")
    print()
  if FLAGS.job_type == "post_treatment":
    print("Parameters for plots:")
    for i in all_params['plots']:
      print(f"{i}: {all_params['plots'][i]}")

  if FLAGS.job_type != 'all':
    all_templates = group_templates(all_params, [FLAGS.job_type])
    create_single_jobfile(FLAGS.job_type, all_templates, run_params, FLAGS.path_to_parameters)
  else:
    all_templates = group_templates(all_params, ['alignment', 'jobarray', 'post_treatment'])
    create_all_jobfile(all_templates, run_params, FLAGS.path_to_parameters)
if __name__ == "__main__":
  app.run(main)
