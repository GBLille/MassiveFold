#!/usr/bin/env python

import json
from string import Template
from absl import app, flags

FLAGS = flags.FLAGS
flags.DEFINE_enum("job_type", 'all', ['all', 'alignment', 'jobarray', 'post_treatment'], 'Type of the jobfile to create')
flags.DEFINE_string("sequence_name", "", 'name of the fasta sequence used for the run.')
flags.DEFINE_string("run_name", '', 'name of the run, it can be anything and it will be the name of the output directory under the sequence name.')
flags.DEFINE_bool("create_files", True, '')
flags.DEFINE_string("path_to_parameters", "", "Path to a json file were the jobfile parameters can be specified.")

def create_single_jobfile(jobfile_type, templates:dict, params):
  jobfile = Template(templates[jobfile_type]).substitute(params)
  if FLAGS.create_files:
    with open(f"{FLAGS.sequence_name}_{FLAGS.run_name}_{jobfile_type}.slurm", 'w') as slurm_job:
      slurm_job.write(jobfile)
  return jobfile

def create_all_jobfile(templates:dict, params:dict):
  for jobtype in ['alignment', 'jobarray', 'post_treatment']:
    create_single_jobfile(jobtype, templates, params)
    
def main(argv):

  with open(f'{FLAGS.sequence_name}_{FLAGS.run_name}_batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
  run_params = {
      'run_name': FLAGS.run_name,
      'sequence_name': FLAGS.sequence_name,
      'substitute_batch_number': list(batches.keys())[-1]
    }

  if not FLAGS.path_to_parameters:
    raise ValueError('Parameters files missing, use --path_to_parameters.')
  with open(FLAGS.path_to_parameters, 'r') as parameters_json:
    all_params = json.load(parameters_json)

  run_params.update(all_params['custom_params'])
  run_params.update(all_params['MF_run'])
  
  if FLAGS.job_type == "jobarray":  
    print("Parameters of the run:")
    for i in all_params['custom_params']:
      print(f"{i}: {all_params['custom_params'][i]}")
    for i in all_params['MF_run']:
      print(f"{i}: {all_params['MF_run'][i]}")

  if 'jeanzay_gpu_memory' in run_params:
    run_params['jeanzay_gpu_memory'] = f"-{run_params['jeanzay_gpu_memory']}"
    run_params['jeanzay_account'] = f"{run_params['jeanzay_project']}@{run_params['jeanzay_gpu']}"
    run_params['jeanzay_full_gpu'] = f"{run_params['jeanzay_gpu']}{run_params['jeanzay_gpu_memory']}"

  grouped_template = all_params['MF_parallel']['grouped_templates']
  with open(grouped_template, 'r') as templates:
    all_templates = json.load(templates)
    
  if FLAGS.job_type != 'all':
    create_single_jobfile(FLAGS.job_type, all_templates, run_params)
  else:
    create_all_jobfile(all_templates, run_params)



if __name__ == "__main__":
  app.run(main)
