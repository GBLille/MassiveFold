#!/usr/bin/env python

from string import Template
import json

if __name__ == "__main__":
  all_templates = {}

  # Store template of each base in a dictionnary
  with open('./templates/alignment_ugsf.slurm', 'r') as file:
    alignment = file.read()
  all_templates['alignment'] = alignment
    
  with open('./templates/jobarray_ugsf.slurm', 'r') as file:
    jobarray = file.read()
  all_templates['jobarray'] = jobarray
  
  with open('./templates/post_treatment_ugsf.slurm', 'r') as file:
    post_treatment = file.read()
  all_templates['post_treatment'] = post_treatment
  
  # Export the templates as a json  
  with open('./templates/templates.json', 'w') as file_out:
    json.dump(all_templates, file_out, indent=4)
