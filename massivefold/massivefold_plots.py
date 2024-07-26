#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import sys
import pickle
import json
from absl import flags
from absl import app
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from plots.colabfold_plots import plot_msa_v2, plot_plddts, plot_confidence, plot_paes, plot_plddt_legend
from scipy.stats import gaussian_kde
import shutil
from matplotlib.lines import Line2D
import seaborn as sns

FLAGS = flags.FLAGS

flags.DEFINE_string('input_path', None, 
                    'Path to directory were the alphafold output to be plotted lies.')
flags.DEFINE_integer('top_n_predictions', 10, 
                     'Specify the number of predictions taken into account for plotting, it will be the n best predictions.')
flags.DEFINE_list('chosen_plots', [], 
                  'Specify the plots you want to get.'
                  'CF_plddt for plddt of each predictions, DM_plddt_PAE for plddt and PAE on the same plot for each prediction,'
                  'CF_PAEs for PAEs of all predictions on the same plot and CF_plddts for plddt of all predictions on the same plot.')
flags.DEFINE_enum('action', "save", ["save", "show"], "Choose to save the plot or show them.")
flags.DEFINE_string('output_path', None, 
                    'Path to directory that will store the plots, same as the input dir by default.')
flags.DEFINE_list('runs_to_compare', [], 'Runs that you want to compare on a same distribution plot')

def extract_top_predictions():
  with open(f'{FLAGS.input_path}/ranking_debug.json', 'r') as json_file:
    scores = json.load(json_file)
  top_pred = scores['order'][:FLAGS.top_n_predictions]
  #top_pred = [pred for pred in top_pred if os.path.isfile(f'{FLAGS.input_path}/result_{pred}.pkl')] 

  return top_pred
  
def CF_PAEs():
  all_models_pae = []
  jobname = FLAGS.input_path
  preds_to_plot = extract_top_predictions()
  
  pkl_dir = jobname
  light_pkl = f'{jobname}/light_pkl'
  if os.path.isdir(light_pkl) and len(os.listdir(light_pkl)) != 0:
    pkl_dir = light_pkl
  
  for pred in preds_to_plot:
    with open(f'{pkl_dir}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    all_models_pae.append(np.asarray(data['predicted_aligned_error']))
  plot_paes(all_models_pae)
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/top_{FLAGS.top_n_predictions}_PAE.png", dpi=200)
    print(f"Saved as top_{FLAGS.top_n_predictions}_PAE.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()
    
def CF_plddts():
  all_models_plddt = []
  jobname = FLAGS.input_path
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    all_models_plddt.append(np.asarray(data['plddt']))
  plot_plddts(all_models_plddt)
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/top_{FLAGS.top_n_predictions}_plddt.png", dpi=200)
    print(f"Saved as top_{FLAGS.top_n_predictions}_plddt.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def MF_DM_dual_plddt_PAE(prediction, rank):
  jobname = FLAGS.input_path
  pkl_dir = jobname
  light_pkl = f'{jobname}/light_pkl'
  if os.path.isdir(light_pkl) and len(os.listdir(light_pkl)) != 0:
    pkl_dir = light_pkl
  with open(f'{pkl_dir}/result_{prediction}.pkl', "rb") as results_file:
    results = pickle.load(results_file)

  pae_outputs = {}
  pae_outputs["test"] = (results["predicted_aligned_error"], results['max_predicted_aligned_error'])
  pae, max_pae = list(pae_outputs.values())[0]
  plt.figure(figsize=[8 * 2, 6])
  
  plt.subplot(1, 2, 1)
  plt.plot(results['plddt'])
  plt.title(f'Predicted LDDT')
  plt.suptitle(f'rank_{rank}_{prediction}')
  plt.xlabel('Residue')
  plt.ylabel('pLDDT')

  plt.subplot(1, 2, 2)
  plt.imshow(pae, vmin=0., vmax=max_pae, cmap='Greens_r')
  plt.colorbar(fraction=0.046, pad=0.04)
  plt.title('Predicted Aligned Error (Ångströms)')
  plt.suptitle(f'rank_{rank}_{prediction}')
  plt.xlabel('Scored residue')
  plt.ylabel('Aligned residue')
  
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/rank_{rank}_{prediction}_plddt_PAE.png", dpi=200)
    print(f"Saved as rank_{rank}_{prediction}_plddt_PAE.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def call_dual():
  preds_to_plot = extract_top_predictions()
  for i, pred in enumerate(preds_to_plot):
    MF_DM_dual_plddt_PAE(pred, i)
  
def MF_indiv_plddt():
  jobname = FLAGS.input_path
  preds_to_plot = extract_top_predictions()
  for i, pred in enumerate(preds_to_plot):
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    plot_confidence(data['plddt'])
    plt.title(f'rank_{i}_{pred} predicted lDDT')
    if FLAGS.action == "save":
      plt.savefig(f"{FLAGS.output_path}/rank_{i}_{pred}_plddt.png", dpi=200)
      print(f"Saved as rank_{i}_{pred}_plddt.png")
      plt.close()
    if FLAGS.action == "show":
      plt.show()
            
def MF_coverage():
  jobname = FLAGS.input_path
  if os.path.isfile(f'{jobname}/features.pkl'):
    with open(f'{jobname}/features.pkl', 'rb') as f:
      data = pickle.load(f)
    plot_msa_v2(data)
    if FLAGS.action == "save":
      plt.savefig(f"{FLAGS.output_path}/alignment_coverage.png", dpi=200)
      print(f"Saved as alignment_coverage.png")
      plt.close()
    elif FLAGS.action == "show":
      plt.show()

def MF_score_histogram(scores:dict):
  try:
    scores = scores['iptm+ptm']
  except KeyError:
    scores = scores['plddts']

  # Global score distribution 
  all_scores = [scores[model] for model in scores]
  histogram, ax1 = plt.subplots()
  ax1.hist(all_scores, bins=50)
  histogram.suptitle('Global score distribution')
  ax1.set(xlabel='Ranking confidence', ylabel='Number of predictions')
  if FLAGS.action == "save":
    histogram.savefig(f"{FLAGS.output_path}/score_distribution.png", dpi=200)
    print("Saved as score_distribution.png")
    plt.close(histogram)

def MF_versions_density(scores:dict):
  try:
    scores = scores['iptm+ptm']
  except KeyError:
    print('\nOnly one version of NN models, no versions density plot.\n')
    return None

  # Score distribution by NN model version
  versions = {prediction.split('multimer_')[1].split('_pred')[0]: [] for prediction in scores}
  
  for model in scores:
    versions[model.split('multimer_')[1].split('_pred')[0]].append(scores[model])
  scores_per_version = pd.DataFrame(versions)
  scores_per_version = scores_per_version.sort_index(axis=1)
  kde_versions, ax2 = plt.subplots()
 
  kde_versions.suptitle('Score density')
  sns.kdeplot(scores_per_version, ax=ax2)
  #scores_per_version.plot(kind="kde", ax=ax2, bw_method=0.3)
  ax2.set(
    xlabel='Ranking confidence',
    ylabel="Density"
  )

  if FLAGS.action == "save":
    kde_versions.savefig(f"{FLAGS.output_path}/versions_density.png",dpi=200)
    print("Saved as versions_density.png")
    plt.close(kde_versions)

def MF_models_scores(scores:dict):
  try:
    scores = scores['iptm+ptm']
    s_type = 'iptm+ptm'
  except KeyError:
    scores = scores['plddts']
    s_type = "plddts"
  # Score distribution by NN model
  NN_models = {prediction.split('_pred')[0]: [] for prediction in scores}
  for model in scores:
    NN_models[model.split('_pred')[0]].append(scores[model])

  scores_per_model = pd.DataFrame(NN_models)
  cols = scores_per_model.columns
  scores_per_model.columns = [col.replace('model_', '').replace('multimer_', '') for col in cols]
  pastel_colors = ['#add8e6', '#87ceeb', '#b0c4de', '#b0e0e6',
  '#e0ffff', '#afeeee', '#90ee90', '#98fb98', '#fafad2', '#ffa07a',
  '#f08080', '#ffb6c1', '#e6a8d7', '#fff0f5', '#ffe4e1']
  colors = {
  'boxes': '#add8e6',
  }
# Create a boxplot with inclined x-axis labels
  fig, ax = plt.subplots()
  ax = scores_per_model.boxplot(sym='g+', patch_artist=True, color = colors, flierprops=dict(markerfacecolor='red'))
  for box, color in zip(ax.artists, pastel_colors):
    box.set_facecolor(color)

  plt.grid(False)
  fig.suptitle("Score per NN model")
  
  ax.set(
    xlabel="NN model",
    ylabel="Ranking confidence"
  )
  ax.set_xticklabels(scores_per_model.columns, rotation=45)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  if s_type == 'iptm+ptm':
    ax.set_ylim(bottom=0, top=1.1)
  elif s_type == 'plddts':
    ax.set_ylim(bottom=0, top=110)
  
  plt.tight_layout()
 
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/models_scores.png", dpi=200)
    print("Saved as models_scores.png")

  if FLAGS.action == "show":
    plt.show()
  plt.close()

def MF_score_distribution():
  distribution_types = ["scores", "versions_scores", "models_scores"]
  jobname = FLAGS.input_path
  with open(f'{jobname}/ranking_debug.json', 'r') as json_scores:
    scores = json.load(json_scores)
  
  DISTRIBUTION_MAP = {
    "scores": MF_score_histogram,
    "versions_scores": MF_versions_density,
    "models_scores":MF_models_scores
  }

  for distrib in distribution_types:
    DISTRIBUTION_MAP[distrib](scores)
    #try:
    #  DISTRIBUTION_MAP[distrib](scores)
    #except:
    #  print(f"Error while trying to plot {DISTRIBUTION_MAP[distrib].__name__}()")

def MF_distribution_comparison():
  sequence_path = FLAGS.input_path
  runs = FLAGS.runs_to_compare
  all_scores = {}
  for run in runs:
    run_basename = os.path.basename(run)
    scores = []
    with open(f"{sequence_path}/{run}/ranking_debug.json" ,'r') as json_file:
      run = json.load(json_file)
      try:
        run_scores = run['iptm+ptm']
      except KeyError:
        run_scores = run['plddts']

    run_scores = [run_scores[score] for score in run_scores]
    if len(runs) <= 3:
      plt.hist(run_scores, bins=30, alpha=0.5, label=run_basename)
    else:
      all_scores[run_basename] = run_scores
  if len(runs) > 3:
    datas = pd.DataFrame(all_scores)
    datas.plot(kind='kde')
  plt.legend()
  plt.title(f'Distribution comparison between {os.path.basename(sequence_path)} runs')
  plt.xlabel('Ranking confidence')
  plt.ylabel('Number of predictions')
  plt.savefig(f'{sequence_path}/distribution_compa.png', dpi=200)
  print('Saved as distribution_compa.png')

def MF_decode_array(encoded_pred:str):
  encoded_lst = [ int(i) for i in encoded_pred.replace('[', '').replace(']', '').split() ]
  if encoded_lst[0]:
    decoded_name = f"model_{encoded_lst[1]}_multimer_v{encoded_lst[2]}_pred_{encoded_lst[3]}"
  else:
    decoded_name = f"model_{encoded_lst[1]}_ptm_pred_{encoded_lst[3]}"

  return decoded_name

def MF_cf_recycling_parser(log_file, prediction, cf_name_map):
  with open(log_file, 'r') as log_recycle:
    lines = log_recycle.readlines()
  with open(cf_name_map, 'r') as json_file:
    name_map = json.load(json_file)
  logs = []
  for line in lines:
    if "--num-recycle" in line:
      max_recycles = int(line.split('--num-recycle')[1].strip())
    if "recycle=" in line:
      recycle_line = line.split(' ')
      for i, elem in enumerate(recycle_line):
        if elem.startswith('alphafold2'):
          recycle_line = recycle_line[i:]
          break
      logs.append(' '.join(recycle_line).strip())

  recycling = {}
  for log in logs:
    elements = log.split(' ')
    prediction = elements[0]
    preset = prediction.split('_')[1]
    pred_name = name_map[prediction]
    if pred_name not in recycling:
      recycling[pred_name] = {
        'scores': [ None for _ in range(max_recycles + 1) ],
        'distances': [ None for _ in range(max_recycles + 1) ]
        }
    n_recycle = int(elements[1].split('recycle=')[1])
    if n_recycle == 0:
       recycling[pred_name]['distances'][0] = 0
    else:
      distance = float((elements[-1].split('tol=')[1]))
      recycling[pred_name]['distances'][n_recycle] = distance

    if preset == 'multimer':
      ptm, iptm = float(elements[3].split('pTM=')[1]), float(elements[4].split('ipTM=')[1])
      recycling[pred_name]['scores'][n_recycle] = 0.8*iptm + 0.2*ptm
    elif preset == 'ptm':
      plddt, ptm = float(elements[2].split('pLDDT=')[1]), float(elements[3].split('pTM=')[1])
      recycling[pred_name]['scores'][n_recycle] = ptm
  return recycling

def MF_recycling_parser(log_file, prediction):
  with open(log_file, 'r') as log_recycle:
    lines = log_recycle.readlines()
  logs = []
  # get recycle nb and its lines from log file
  for line in lines:
    if '--max_recycles=' in line:
      max_recycles = int(line.split('--max_recycles=')[1].strip())
    if line.startswith('['):
      logs.append(line.strip())

  encoded_recycles = {}
  for log in logs:
    if 'recycle' in log:
      n_recycle = int(log.split('recycle=')[1].split(' ')[0])
    elif 'last step' in log:
      n_recycle = 'last'
    else:
      print('Unexpected log')
      print(log)
      print('Pass this log')
      continue

    pred_encoded = log.split('] ')[0][1:]
    pred_name = MF_decode_array(pred_encoded)
    if pred_name not in encoded_recycles:
      encoded_recycles[pred_name] = {
        'scores': [ None for _ in range(max_recycles) ],
        'distances': [ None for _ in range(max_recycles) ],
        'last': {}
        }
    
    if 'distance' in log:
      distance = log.split('distance=')[1].split(' ')[0]
      if n_recycle == 0:
        encoded_recycles[pred_name]['distances'][0] = 0
      elif n_recycle != 1:
        try:
          encoded_recycles[pred_name]['distances'][n_recycle - 1] = float(distance)
        except TypeError:
          encoded_recycles[pred_name]['last']['distance'] = float(distance)
    elif 'confidence' in log:
      score = log.split('confidence=')[1].split(' ')[0]
      try:
        encoded_recycles[pred_name]['scores'][n_recycle] = float(score)
      except TypeError:
        encoded_recycles[pred_name]['last']['score'] = float(score)

  for pred in encoded_recycles:
    single_pred = encoded_recycles[pred]
    lasts = single_pred['last']
    
    none_index = single_pred['scores'].index(None) if None in single_pred['scores'] else len(single_pred['scores'])
    single_pred['scores'].insert(none_index, lasts['score'])

    none_index = single_pred['distances'].index(None) if None in single_pred['distances'] else len(single_pred['distances'])
    single_pred['distances'].insert(none_index, lasts['distance'])
    del single_pred['last']

  return encoded_recycles

def MF_recycling_export(all_predictions, file):
  with open(file, 'w') as json_output:
    json.dump(all_predictions, json_output, indent=4)

def MF_recycling_pred_to_batch(batches_file):
  with open(batches_file, 'r') as batches:
    batches_to_model = json.load(batches)
  model_to_batch = {}
  for batch in batches_to_model:
    start = int(batches_to_model[batch]['start'])
    end = int(batches_to_model[batch]['end'])
    model = batches_to_model[batch]['model']
    if not model in model_to_batch:
      model_to_batch[model] = {}
    pred_to_batch = {pred: batch for pred in range(start, end + 1)}
    model_to_batch[model].update(pred_to_batch)
  return model_to_batch

def MF_recycling_plot(
    seq:str,
    run:str,
    all_values:dict,
    prediction_to_plot:str,
    rank:int, 
    tol:float,
    preset:str):
 
  recycle_dir = f'{FLAGS.output_path}/recycles/'
  if not shutil.os.path.exists(recycle_dir):
    shutil.os.makedirs(recycle_dir)

  fig, ax1 = plt.subplots()
  scores_elements = all_values[prediction_to_plot]['scores']
  distances_elements = all_values[prediction_to_plot]['distances']
  n_recycling = np.linspace(0, len(scores_elements) - 1, len(scores_elements))

  color = 'tab:grey'
  ax1.plot(n_recycling, [tol for val in n_recycling], color=color)
  ax1.set_xlabel('Recycling steps')
  tol_pos_increment = max([ i for i in distances_elements if i ])/50
  plt.text(n_recycling[0], tol + tol_pos_increment, 'Early stop tolerance', color = color)
 
  ax2 = ax1.twinx()
  color = 'tab:red'
  if preset == 'multimer':
    metric = 'Ranking confidence'
  elif preset == 'ptm':
    metric = 'pTM'
  ax2.set_ylabel(metric, color=color)
  ax2.set_ylim(bottom=0, top=1.1)
  ax2.plot(n_recycling, scores_elements, color=color, alpha=0.8)
  ax2.tick_params(axis='y', labelcolor=color)

  color = 'tab:blue'
  ax1.set_ylabel('Distance with previous step structure', color=color)
  ax1.plot(n_recycling, distances_elements, color=color, alpha=0.8)
  ax1.tick_params(axis='y', labelcolor=color)

  locator = ticker.MaxNLocator(integer=True)
  ax1.xaxis.set_major_locator(locator)
  ax2.xaxis.set_major_locator(locator)

  fig.suptitle(f"{seq} - {run}")
  plt.title(f"ranked_{rank}_{prediction_to_plot}")
  fig.tight_layout()

  if FLAGS.action == "save":
    plt.savefig(f"{recycle_dir}/ranked_{rank}_unrelaxed_{prediction_to_plot}.png", dpi=200)
    print(f"Saved as recycles/ranked_{rank}_unrelaxed_{prediction_to_plot}.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def MF_recycling():
  all_models_pae = []
  jobname = FLAGS.input_path 
  run, seq = os.path.basename(jobname), os.path.basename(os.path.dirname(jobname))
  logs_dir = os.path.join(jobname, '../../../log/', seq, run)
  batches_file = os.path.join(logs_dir, f'{seq}_{run}_batches.json')
  
  preds_to_plot = extract_top_predictions()
  model_to_batch = MF_recycling_pred_to_batch(batches_file) 
  
  all_recycling_values = {}
  for pred in preds_to_plot:
    if pred not in all_recycling_values:
      pred_model, pred_nb = pred.split('_pred_')
      model_preds = model_to_batch[pred_model]
      batch_number = model_preds[int(pred_nb)]
      log_file = os.path.join(logs_dir, f"jobarray_{batch_number}.log")
      # parse logs depending on inference engine specificity
      colabfold_name_map = os.path.join(jobname, 'unified_map.json')
      if os.path.isfile(colabfold_name_map):
        pred_recycling = MF_cf_recycling_parser(log_file, pred, colabfold_name_map)
      else:
        pred_recycling = MF_recycling_parser(log_file, pred)
      all_recycling_values.update(pred_recycling)

  with open(os.path.join(logs_dir, 'jobarray_0.log'), 'r') as logfile:
    lines = logfile.readlines()

  # parse tol and preset depending on AFmassive or ColabFold
  try:
    early_stop_tolerance = float([line.split('=')[1].strip() for line in lines if 'early_stop_tolerance=' in line][0])
    preset = 'multimer'
  except IndexError:
    early_stop_tolerance= float([line.split('tolerance')[1].strip() for line in lines if '--recycle-early-stop-tolerance' in line][0])
    preset = [ line.split('alphafold2_')[1].split('_')[0] for line in lines if 'alphafold2_' in line ][0].strip()

  recycling_file = os.path.join(jobname, 'recycling_log.json')
  MF_recycling_export(all_recycling_values, recycling_file)
  with open(recycling_file, 'r') as recycling_json:
    recycling = json.load(recycling_json)
  
  for i, pred in enumerate(preds_to_plot):
    MF_recycling_plot(
      seq,
      run,
      recycling,
      pred,
      i,
      early_stop_tolerance,
      preset)

def main(argv):
  FLAGS.input_path = os.path.realpath(FLAGS.input_path)
  MF_plots = {
    "DM_plddt_PAE": call_dual,
    "CF_plddt": MF_indiv_plddt,
    "CF_PAEs": CF_PAEs,
    "CF_plddts": CF_plddts,
    "coverage": MF_coverage,
    "score_distribution": MF_score_distribution,
    "distribution_comparison": MF_distribution_comparison,
    "recycles": MF_recycling
    }

  # Flags checking
  if not FLAGS.input_path or not FLAGS.chosen_plots:
    print('Required flags: --input_path and --chosen_plots')  
  if not FLAGS.output_path:
    FLAGS.output_path = f"{FLAGS.input_path}/plots/"
  print(f"Plot are stored in {FLAGS.output_path}")
  if "distribution_comparison" in FLAGS.chosen_plots and not FLAGS.runs_to_compare:
    print('Flag --runs_to_compare is required for --chosen_plots=distribution_comparison')
    FLAGS.chosen_plots = [ plot for plot in FLAGS.chosen_plots if plot != 'distribution_comparison' ]

  if not shutil.os.path.exists(FLAGS.output_path) and 'distribution_comparison' not in FLAGS.chosen_plots:
    shutil.os.makedirs(FLAGS.output_path)
  for chosen_plot in FLAGS.chosen_plots:
    MF_plots[chosen_plot]()

if __name__ == "__main__":
  """ 
  These functions are a combination of MassiveFold's team work and the following scripts from ColabFold repository and DeepMind colab notebook:
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py
  https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb
  
  Here are some basic commands:
  python MF_plots.py --input_path ./jobname --chosen_plots coverage,CF_PAEs
    -> regardless of the plot type, plot alignment coverage and group PAE for top 10 predictions
  """
  app.run(main)
