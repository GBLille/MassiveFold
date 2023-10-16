#!/usr/bin/env python

import math
import json
import numpy as np
from copy import deepcopy
from absl import app, flags



flags.DEFINE_integer('predictions_per_model', 25, 
                     'Choose the number of predictions inferred by each neural network model.')
flags.DEFINE_integer('batch_size', 25, 
                     'Standard size of a prediction batch, if the number of prediction per model\
                       is not a multiple of it, the last batch will be smaller .')
flags.DEFINE_list('models_to_use', None, 'Select the models used for prediction among the five models of each AlphaFold2 version (15 in total).')

FLAGS = flags.FLAGS

def batches_per_model(pred_nb_per_model:int):
  opt_batch_nb = math.ceil(pred_nb_per_model/FLAGS.batch_size)
  batch_sizes = []
  for _ in range(1, opt_batch_nb+1):
    # split total by batch of the same size, the remaining is a single smaller batch
    if pred_nb_per_model - FLAGS.batch_size >= 0:
      pred_nb_per_model -= FLAGS.batch_size
      batch_sizes.append(FLAGS.batch_size)
    else:
      #print(f"Last batch size: {pred_nb_per_model}")
      batch_sizes.append(pred_nb_per_model)
  
  batch_edges = list(np.cumsum([0] + batch_sizes))
  one_model_batches = {i+1: {'start': str(batch_edges[i]), 'end': str(batch_edges[i+1]-1)} for i in range(len(batch_edges)-1)}
  return one_model_batches

def batches_all_models(batches_unit, all_models):
  batches = {}

  for i, model in enumerate(all_models):
    unadded_batch = deepcopy(batches_unit)
    for batch in unadded_batch:
      unadded_batch[batch].update({'model': model})
      batches[str(batch+i*len(batches_unit)-1)] = unadded_batch[batch]
      
  return batches
      
def main(argv):
  model_names = [
  'model_1_multimer_v1',
  'model_2_multimer_v1',
  'model_3_multimer_v1',
  'model_4_multimer_v1',
  'model_5_multimer_v1',
  'model_1_multimer_v2',
  'model_2_multimer_v2',
  'model_3_multimer_v2',
  'model_4_multimer_v2',
  'model_5_multimer_v2',
  'model_1_multimer_v3',
  'model_2_multimer_v3',
  'model_3_multimer_v3',
  'model_4_multimer_v3',
  'model_5_multimer_v3'
  ]

  FLAGS.models_to_use = [model for model in model_names if model in FLAGS.models_to_use]
  if not FLAGS.models_to_use:
    FLAGS.models_to_use = model_names
    
  per_model_batches = batches_per_model(pred_nb_per_model=FLAGS.predictions_per_model)
  all_model_batches = batches_all_models(per_model_batches, FLAGS.models_to_use)
  
  with open("batches.json", "w") as json_output:
    json.dump(all_model_batches, json_output)
      
if __name__ == "__main__":
  app.run(main)
