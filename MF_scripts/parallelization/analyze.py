#!/usr/bin/env python

from absl import flags, app
import pandas as pd
import json
import os 

FLAGS = flags.FLAGS
flags.DEFINE_string("input", "", "")
flags.DEFINE_integer("top_n", 5, "")
flags.DEFINE_enum("get", "models", ["models", "analytics"], "")


def order_mean(data2):
  ordered = sorted(dict(data2).items(), key=lambda x:x[1], reverse=True)
  df = {"object": [], "mean": []}
  for model in ordered:
    df['object'].append(model[0])
    df['mean'].append(model[1])
  print(pd.DataFrame(df))

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
    if model  not in models:
      models[model] = scores['iptm+ptm'][prediction]
    if len(models) == number:
      break
  print(','.join(list(models)))

def main(argv):
  path = os.path.expanduser(FLAGS.input)
  
  with open(f"{path}/ranking_debug.json", "r") as data:
    data = json.load(data)
  
  if FLAGS.get == "models":
    get_top_models(data, FLAGS.top_n)
  elif FLAGS.get == "analytics":
    some_stats(data)


if __name__ == "__main__":
  app.run(main)
