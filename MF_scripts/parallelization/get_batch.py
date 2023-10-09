#!/usr/bin/env python

from absl import app, flags
import json

FLAGS = flags.FLAGS

flags.DEFINE_string('batch_id', '', 'Task id or batch number from which the element are retrieved.')
flags.DEFINE_string('element', '',
                  'Select the element (start, end or model) that is retrieved by the script.')
flags.DEFINE_string('json_path', './batches.json', 'Path of the json from which the element are retrieved.')

def main(argv):
  if not FLAGS.batch_id or not FLAGS.element:
    raise ValueError('\nUsage: ./get_batch.py --batch_id [str:batch_id] --element [str:element] [[options]]')
  if FLAGS.element not in ['start', 'end', 'model']:
    raise ValueError('--element is either start, end or model.')
  
  with open(FLAGS.json_path, 'r') as f:
    all_batches = json.load(f)
    
  try:
    print(all_batches[FLAGS.batch_id][FLAGS.element])
  except KeyError:
    raise ValueError(f"Either no batch {FLAGS.batch_id} or no element '{FLAGS.element}' for this batch.")

if __name__ == "__main__":
  app.run(main)