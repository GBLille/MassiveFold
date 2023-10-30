#!/usr/bin/env python

from absl import flags, app
from string import Template
import json

FLAGS = flags.FLAGS
flags.DEFINE_string('parameters', '', 'Parameters file containing the template paths')

def group_template(all_template_paths):
  grouped_templates = {}                                                               
                                                                                   
  # Store template of each base in a dictionnary                                   
  with open(all_template_paths['alignment'], 'r') as file:     
    alignment = file.read()                                                        
  grouped_templates['alignment'] = alignment                                           
                                                                                   
  with open(all_template_paths['jobarray'], 'r') as file:      
    jobarray = file.read()                                                         
  grouped_templates['jobarray'] = jobarray                                             
                                                                                   
  with open(all_template_paths['post_treatment'], 'r') as file:
    post_treatment = file.read()                                                   
  grouped_templates['post_treatment'] = post_treatment                                 
                                                                                   
  # Export the templates as a json                                                 
  with open(all_template_paths['group'], 'w') as file_out:  
    json.dump(grouped_templates, file_out, indent=4)                                   


def main(argv):
  with open(FLAGS.parameters, 'r') as params:
    MF_run_params = json.load(params)['MF_parallel']
  
  template_paths = {
    "alignment": MF_run_params['alignment_template'],
    "jobarray": MF_run_params['jobarray_template'],
    "post_treatment": MF_run_params['post_treatment_template'],
    "group": MF_run_params['grouped_templates']
}

  group_template(template_paths)

if __name__ == "__main__":
  app.run(main)
