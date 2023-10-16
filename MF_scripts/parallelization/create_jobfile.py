#!/usr/bin/env python

import json
from string import Template
from absl import app, flags

FLAGS = flags.FLAGS
flags.DEFINE_enum("job_type", 'all', ['all', 'alignment', 'jobarray', 'post_treatment'], 'Type of the jobfile to create')
flags.DEFINE_enum("cluster_name", 'ugsf', ['ugsf', 'jeanzay'], 'Specify the name of the cluster where the job is run.')
flags.DEFINE_string("jobname", '', 'name of the job and fasta file as the input sequence')
flags.DEFINE_bool("create_files", True, '')

def create_single_jobfile(jobfile_type, templates:dict, params):
  jobfile = Template(templates[jobfile_type]).substitute(params)
  if FLAGS.create_files:
    with open(f"{FLAGS.jobname}_{jobfile_type}.slurm", 'w') as slurm_job:
      slurm_job.write(jobfile)
  return jobfile

def create_all_jobfile(templates:dict, params:dict):
  for jobtype in ['alignment', 'jobarray', 'post_treatment']:
    create_single_jobfile(jobtype, templates, params)
    
def main(argv):
  params = {}

  with open('batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
    params.update({
      'jobname': FLAGS.jobname,
      'substitute_batch_number': list(batches.keys())[-1]
    })

  if FLAGS.cluster_name == 'jeanzay':
    # to do: replace this by a json option file
    custom_params = {               
      'jeanzay_gpu': 'v100',      
      'jeanzay_gpu_memory': '16g',
      'jeanzay_project': 'fvp',   
      }                      
    params.update(custom_params)

    params['jeanzay_account'] = f"{params['jeanzay_project']}@{params['jeanzay_gpu']}"
    params['jeanzay_full_gpu'] = f"{params['jeanzay_gpu']}-{params['jeanzay_gpu_memory']}"
  
  with open('batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
    params['substitute_batch_number'] = list(batches.keys())[-1]
  
  with open(f'./templates/templates_{FLAGS.cluster_name}.json', 'r') as templates:
    all_templates = json.load(templates)
    
  if FLAGS.job_type != 'all':
    create_single_jobfile(FLAGS.job_type, all_templates, params)
  else:
    create_all_jobfile(all_templates, params)

if __name__ == "__main__":
  app.run(main)
