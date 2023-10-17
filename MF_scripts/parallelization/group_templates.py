#!/usr/bin/env python

from absl import flags, app
from string import Template
import json

FLAGS = flags.FLAGS
flags.DEFINE_enum('cluster_name', 'all', ['all', 'ugsf', 'jeanzay', 'generic'], 'The cluster on which the job is run.')

def group_template(cluster_name):
  all_templates = {}                                                               
                                                                                   
  # Store template of each base in a dictionnary                                   
  with open(f'./templates/alignment_{cluster_name}.slurm', 'r') as file:     
    alignment = file.read()                                                        
  all_templates['alignment'] = alignment                                           
                                                                                   
  with open(f'./templates/jobarray_{cluster_name}.slurm', 'r') as file:      
    jobarray = file.read()                                                         
  all_templates['jobarray'] = jobarray                                             
                                                                                   
  with open(f'./templates/post_treatment_{cluster_name}.slurm', 'r') as file:
    post_treatment = file.read()                                                   
  all_templates['post_treatment'] = post_treatment                                 
                                                                                   
  # Export the templates as a json                                                 
  with open(f'./templates/templates_{cluster_name}.json', 'w') as file_out:  
    json.dump(all_templates, file_out, indent=4)                                   



def main(argv):
  pipe_separated_enum_list = flags.FLAGS['cluster_name'].help.split(':')[0][1:-1]
  clusters_enum_list = pipe_separated_enum_list.split('|')
 
  if FLAGS.cluster_name != 'all':
    group_template(FLAGS.cluster_name)
  else:
    clusters_enum_list.remove('all')
    for cluster in clusters_enum_list:
      group_template(cluster)

if __name__ == "__main__":
  app.run(main)
