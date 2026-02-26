#!/usr/bin/env python

from absl import flags, app
import pandas as pd
from math import floor, ceil
import json
import os 

FLAGS = flags.FLAGS
flags.DEFINE_enum(
  "get", "models", ["models", "batch_calibration", "analytics"],
  "Specify the role of the script."
  )
flags.DEFINE_string(
  "input", "",
  "When --get='models', --input designates the output directory containing"
  "the ranking file where --top_n models are extracted from"
  "When --get='batch_calibration', --input designates the log directory path"
  "of the run containing the times taken for prediction"
  )
flags.DEFINE_integer(
  "top_n", 5, "The number of top AF2 NN models to be extracted"
  "Specify how many of the best models should be included in the final selection."
  )
flags.DEFINE_float('wall_time', 20, 'Inference time in hour to not exceed.')
flags.DEFINE_float('add_excess', 0.1, 'Excess time proportion for the inference of a single prediction')

def extract_longer(jobarray_path):
  time_lines = os.popen(f"cat {jobarray_path}/jobarray_* | grep 'predict time'").read()
  if not time_lines:
    print(f'No prediction for {os.path.basename(FLAGS.logs_dir)} yet')
    exit()
  i_thing = time_lines.split(' ')[0]
  time_lines_list = time_lines.split(i_thing)

  lines = [line[1:-1] for line in time_lines_list]
  times = [float(line.split(' ')[-1].replace('s', '')) for line in lines[1:]]
  return ceil(max(times))

def max_pred_nb_for_walltime(wall_time, single_pred_time):
  wall_time_m = wall_time*60
  wall_time_s = wall_time_m*60
  max_number = wall_time_s/single_pred_time
  return floor(max_number)

def safen_time(single_pred_time, excess_proportion):
  return floor(single_pred_time + (single_pred_time*excess_proportion))

def some_stats(data):
  data_dict = {"prediction": [], "order": [], "score": []}
  for order, prediction in enumerate(data['order']):
    data_dict["prediction"].append(prediction)
    data_dict['order'].append(order)
    data_dict['score'].append(data['iptm+ptm'][prediction])
  
  df = pd.DataFrame(data_dict)
  df['model'] = df['prediction'].str.split('_pred').str[0]
  df['version'] = df['prediction'].str.slice(start=17, stop=19)
  
  model_means = df.groupby('model')['score'].mean().reset_index()
  model_means.columns = ['models', 'scores']
  model_means = model_means.sort_values(by='scores', ascending=False)
  print(model_means)

  version_means = df.groupby('version')['score'].mean()
  
def get_top_models(scores, number):
  models = {}
  # top models selection
  for prediction in scores['order']:
    model = prediction.split('_pred')[0]
    if model not in models:
      models[model] = scores['iptm+ptm'][prediction]
    if len(models) == number:
      break
  print(','.join(list(models)))

def main(argv):
  if FLAGS.get == "models":
    path = os.path.expanduser(FLAGS.input)
    with open(f"{path}/ranking_debug.json", "r") as data:
      data = json.load(data)  
    get_top_models(data, FLAGS.top_n)

  elif FLAGS.get == "batch_calibration":
    path_to_logs = FLAGS.input
    wall_time_h = FLAGS.wall_time
    max_time = extract_longer(path_to_logs)
    safe_time = safen_time(max_time, FLAGS.add_excess)
    print(max_pred_nb_for_walltime(wall_time_h, safe_time))

  elif FLAGS.get == "analytics":
    some_stats(data)

if __name__ == "__main__":
  app.run(main)
