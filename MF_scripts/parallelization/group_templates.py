#!/usr/bin/env python

from absl import flags, app
from string import Template
import json

FLAGS = flags.FLAGS
flags.DEFINE_enum('cluster_name', 'ugsf', ['ugsf', 'jeanzay'], 'The cluster on which the job is run.')

def main(argv):
  all_templates = {}

  # Store template of each base in a dictionnary
  with open(f'./templates/alignment_{FLAGS.cluster_name}.slurm', 'r') as file:
    alignment = file.read()
  all_templates['alignment'] = alignment
    
  with open(f'./templates/jobarray_{FLAGS.cluster_name}.slurm', 'r') as file:
    jobarray = file.read()
  all_templates['jobarray'] = jobarray
   
  with open(f'./templates/post_treatment_{FLAGS.cluster_name}.slurm', 'r') as file:
    post_treatment = file.read()
  all_templates['post_treatment'] = post_treatment 

  # Export the templates as a json  
  with open(f'./templates/templates_{FLAGS.cluster_name}.json', 'w') as file_out:
    json.dump(all_templates, file_out, indent=4)


if __name__ == "__main__":
  app.run(main)
