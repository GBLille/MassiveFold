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
from colabfold_plots import plot_msa_v2, plot_plddts, plot_confidence, plot_paes, plot_plddt_legend
from scipy.stats import gaussian_kde
import shutil


FLAGS = flags.FLAGS

flags.DEFINE_string('input_path', None, 
                    'Path to directory were the alphafold output to be plotted lies.')
flags.DEFINE_integer('top_n_predictions', 10, 
                     'Specify the number of predictions taken into account for plotting, it will be the n best predictions.')
flags.DEFINE_list('chosen_plots', [], 
                  'Specify the plots you want to get.'
                  'CF_plddt for plddt of each predictions, DM_plddt_PAE for plddt and PAE on the same plot for each prediction,'
                  'CF_PAEs for all predictions PAEs on the plot and CF_plddts for all predictions plddt on the same plot.')
flags.DEFINE_enum('action', "save", ["save", "show"], "Chose to save the plot or show them.")
flags.DEFINE_string('output_path', None, 
                    'Path to directory that will store the plot, same as the input dir by default.')
flags.DEFINE_list('runs_to_compare', [], 'Runs you want to compare on their distribution')

def extract_top_predictions():
  with open(f'{FLAGS.input_path}/ranking_debug.json', 'r') as json_file:
    scores = json.load(json_file)
  top_pred = scores['order'][:FLAGS.top_n_predictions]
  top_pred = [pred for pred in top_pred if os.path.isfile(f'{FLAGS.input_path}/result_{pred}.pkl')]
  return top_pred
  
def CF_PAEs():
  all_models_pae = []
  jobname = FLAGS.input_path
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    all_models_pae.append(np.asarray(data['predicted_aligned_error']))
  plot_paes(all_models_pae)
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/top_{FLAGS.top_n_predictions}_PAE.png")
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
    plt.savefig(f"{FLAGS.output_path}/top_{FLAGS.top_n_predictions}_plddt.png")
    print(f"Saved as top_{FLAGS.top_n_predictions}_plddt.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def MF_DM_dual_plddt_PAE(prediction, rank):
  jobname = FLAGS.input_path
  with open(f'{jobname}/result_{prediction}.pkl', "rb") as results_file:
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
    plt.savefig(f"{FLAGS.output_path}/rank_{rank}_{prediction}_plddt_PAE.png")
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
      plt.savefig(f"{FLAGS.output_path}/rank_{i}_{pred}_plddt.png")
      print(f"Saved as rank_{i}_{pred}_plddt.png")
      plt.close()
    if FLAGS.action == "show":
      plt.show()
            
def MF_coverage():
  jobname = FLAGS.input_path
  with open(f'{jobname}/features.pkl', 'rb') as f:
    data = pickle.load(f)
  plot_msa_v2(data)
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/alignment_coverage.png")
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
  ax1.set(title=f"Histogram of {os.path.basename(FLAGS.input_path)}'s \
{len(all_scores)} predictions score distribution", xlabel='Ranking confidence',
ylabel='Prediction number')
  if FLAGS.action == "save":
    histogram.savefig(f"{FLAGS.output_path}/score_distribution.png")
    print("Saved as score_distribution.png")
    plt.close(histogram)

def MF_versions_density(scores:dict):
  try:
    scores = scores['iptm+ptm']
  except KeyError:
    print('\nOnly one version of NN models, no versions density plot.\n')
    return None

  # Score distribution by NN model version
  available_version = {model.split('multimer_')[1].split('_pred')[0] for model in scores}
  scores_per_version = pd.DataFrame(
    {
    'v1': [scores[model] for model in scores if "v1" in model],
    'v2': [scores[model] for model in scores if "v2" in model],
    'v3': [scores[model] for model in scores if "v3" in model],
    }
  )
  kde_versions, ax2 = plt.subplots()
  scores_per_version.plot(kind="kde", ax=ax2, bw_method=0.3)
  ax2.set(
    title="Ranked confidence density per NN model version",
    xlabel="Ranked confidence",
    ylabel="Density"
  )
  if FLAGS.action == "save":
    kde_versions.savefig(f"{FLAGS.output_path}/versions_density.png")
    print("Saved as versions_density.png")
    plt.close(kde_versions)

def MF_models_scores(scores:dict):
  try:
    scores = scores['iptm+ptm']
  except KeyError:
    scores = scores['plddts']

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
  
  ax.set(
    title="Ranking confidence boxplot per NN model",
    xlabel="NN model",
    ylabel="Ranking confidence"
  )
  ax.set_xticklabels(scores_per_model.columns, rotation=45)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.set_ylim(bottom=0, top=1.1)
  plt.tight_layout()
  
  if FLAGS.action == "save":
    plt.savefig(f"{FLAGS.output_path}/models_scores.png", dpi=100)
    print("Saved as models_scores.png")

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

def MF_distribution_comparison():
  sequence_path = FLAGS.input_path
  runs = FLAGS.runs_to_compare
  all_scores = {}
  for run in runs:
    run_basename = os.path.basename(run)
    scores = []
    with open(f"{sequence_path}/{run}/ranking_debug.json" ,'r') as json_file:
      run_scores = json.load(json_file)['iptm+ptm']
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
  plt.savefig(f'{sequence_path}/distribution_compa.png')
  print('Saved as distribution_compa.png')

def MF_extract_pred_recycle(log, relative_position):
  with open(log, 'r') as log_file:
    lines = log_file.readlines()
  recycle_scores = [ float(score.strip()) for score in lines if score.strip()[:2] in ['0', '1', '0.'] ]
  start_symbol = 0
  end_symbol = 1

  validity_conditions = [
    recycle_scores[0] == start_symbol,
    recycle_scores[-1] == end_symbol,
    all([ recycle_scores[index+1] == start_symbol for index, element in enumerate(recycle_scores[:-1]) if element == end_symbol ])
  ]
  
  if all(validity_conditions):
    splited_recycle_scores = []
    i = 0
    for score in recycle_scores[:]:
      if score == 1:
        splited_recycle_scores.append(recycle_scores[1:i])
        recycle_scores = recycle_scores[i+1:]
        i = 0
      else:
        i+=1
    return splited_recycle_scores[relative_position[0]]

def MF_recycles():
  all_models_pae = []
  jobname = FLAGS.input_path
  run = os.path.basename(jobname)
  seq = os.path.basename(os.path.dirname(jobname))
  logs_dir = os.path.join(jobname, '../../../log_parallel/', seq, run)
  batches_file = os.path.join(logs_dir, f'{seq}_{run}_batches.json')
  
  preds_to_plot = extract_top_predictions()

  recycle_dir = f'{FLAGS.output_path}/recycles/'
  if not shutil.os.path.exists(recycle_dir):
    shutil.os.makedirs(recycle_dir)
    
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
  
  for pred in preds_to_plot:
    pred_model, pred_nb = pred.split('_pred_')
    
    model_preds = model_to_batch[pred_model]
    batch_number = model_preds[int(pred_nb)]
    log_file = os.path.join(logs_dir, f"jobarray_{batch_number}.log")
    all_preds_in_batch = sorted([ pred for pred in model_preds if model_preds[pred] == batch_number ])
    position_in_batch = all_preds_in_batch.index(int(pred_nb))
    batch_size = len(all_preds_in_batch)

    relative_position = (position_in_batch, batch_size)

    pred_recycles = MF_extract_pred_recycle(log_file, relative_position)
    
    if pred_recycles:
      fig, ax = plt.subplots()
    
      ax.plot(np.linspace(0, len(pred_recycles) +1, len(pred_recycles)), pred_recycles)
      ax.set_title(f"{pred_model}_pred_{pred_nb} score evolution during recycles")
      ax.set_ylabel('Ranking confidence')
      ax.set_xlabel('Recycles')
      ax.set_ylim(bottom=0, top=1.1)

      if FLAGS.action == "save":
        plt.savefig(f"{recycle_dir}/{pred_model}_pred_{pred_nb}.png")
        print(f"Saved as {pred_model}_pred_{pred_nb}.png")
        plt.close()
      if FLAGS.action == "show":
        plt.show()
    else:
      print(f'Recycles are broken for {os.path.basename(log_file)}')

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
    "recycles": MF_recycles
    }

  # Flags checking
  if not FLAGS.input_path or not FLAGS.chosen_plots:
    print('Required flags: --input_path and --chosen_plots')  
  if not FLAGS.output_path:
    FLAGS.output_path = f"{FLAGS.input_path}/plots/"
  if "distribution_comparison" in FLAGS.chosen_plots and not FLAGS.runs_to_compare:
    print('Flag --runs_to_compare is required for --chosen_plots=distribution_comparison')
    FLAGS.chosen_plots = [ plot for plot in FLAGS.chosen_plots if plot != 'distribution_comparison' ]


  if not shutil.os.path.exists(FLAGS.output_path):
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
