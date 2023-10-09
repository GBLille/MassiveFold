#!/usr/bin/env python

import json
from string import Template

if __name__ == "__main__":
  with open('batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
  
  params = {
    'substitute_batch_number': list(batches.keys())[-1]
  }
  with open('./template_ugsf.slurm', 'r') as f:
    src = Template(f.read())
    result = src.substitute(params)
    print(result)
