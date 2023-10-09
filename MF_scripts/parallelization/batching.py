#!/usr/bin/env python

import math
import json
import numpy as np
from copy import deepcopy
from absl import app, flags


FLAGS = flags.FLAGS

flags.DEFINE_integer('predictions_per_model', 25, 
                     'Choose the number of predictions inferred by each neural network model.')
flags.DEFINE_integer('batch_size', 25, 
                     'Standard size of a prediction batch, if the number of prediction per model\
                       is not a multiple of it, the last batch will be smaller .')
flags.DEFINE_list('models_to_use', ['model_1_multimer_v1',
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
                                    'model_5_multimer_v3'], 
                    'Select the models used for prediction among the five models of each AlphaFold2 version (15 in total).')
flags.DEFINE_boolean('save_json', True, 'Save the batches resulting from the split in a json file named batches.json.')

def batches_per_model(pred_nb_per_model:int):
  opt_batch_nb = math.ceil(pred_nb_per_model/FLAGS.batch_size)
  #print(f"Number of batch per model: {opt_batch_nb}")
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
      
def create_jobfiles():
  pass

def main(argv):
  if not set(FLAGS.models_to_use).issubset(flags.FLAGS['models_to_use'].default):
    raise ValueError(f"\n--models_to_use arguments must be choosed from this list:\
      \n{flags.FLAGS['models_to_use'].default}\n\nIf you want to use all previous models don't specify the --models_to_use option.")

  per_model_batches = batches_per_model(pred_nb_per_model=FLAGS.predictions_per_model)
  all_model_batches = batches_all_models(per_model_batches, FLAGS.models_to_use)
  
  if FLAGS.save_json:
    with open("batches.json", "w") as json_output:
      json.dump(all_model_batches, json_output)
  else:
    for i in all_model_batches:
      print(f"Batch {i}: {all_model_batches[i]}")
      
if __name__ == "__main__":
  app.run(main)