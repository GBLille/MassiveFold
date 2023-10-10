#!/usr/bin/env python

import json
from string import Template
from absl import app, flags

FLAGS = flags.FLAGS
flags.DEFINE_string("jobname", '', 'name of the job and fasta file as the input sequence')
flags.DEFINE_boolean("with_output", False, 'create the jobfile to organize the output and create plot from it.')

def main(argv):
  with open('batches.json', 'r') as json_batches:
    batches = json.load(json_batches)
  params = {
    'substitute_batch_number': list(batches.keys())[-1]
    }
    
  with open('./template_ugsf_jobarray.slurm', 'r') as f:
    src = Template(f.read())
    result = src.substitute(params)
    
  if not FLAGS.jobname:
    print(result)
  else:
    with open(f"{FLAGS.jobname}.slurm", 'w') as jobfile:
      jobfile.write(result)
      
if __name__ == "__main__":
  app.run(main)
