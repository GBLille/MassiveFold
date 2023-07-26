import numpy as np
import os
import sys
import pickle
import json
from absl import flags
from absl import app
import matplotlib.pyplot as plt
from colabplots import plot_msa_v2, plot_plddts, plot_confidence, plot_paes, plot_plddt_legend

FLAGS = flags.FLAGS

flags.DEFINE_string('jobname', None, 'Path to directory were the alphafold output to be plotted lies.')
flags.DEFINE_enum('plot_type', "None", ['one_for_all', 'for_each', "None"],
                  'Select either a global plot for the specified prediction or one plot for each of them.')
flags.DEFINE_integer('n_top_predictions', 10, 'Specify the number of prediction taken into account in the plots, it will be the n best predictions.')
flags.DEFINE_list('chosen_plots', [], 'Specify the plots you want to get.')
flags.DEFINE_enum('action', "save", ["save", "show"], "Chose to save the plot or show them.")
flags.DEFINE_string('output', None, 'Path to directory that will be store the plot.')

def plot_for_each():
  if FLAGS.chosen_plots:
    chosen_types = [FLAGS.chosen_plots]
  else:
    chosen_types = ["DM_dual", "CF_plddt"]
  type_for_each = {"DM_dual": call_dual,
                   "CF_plddt": MF_indiv_plddt}
  for plot in chosen_types:
    type_for_each[plot]()

def extract_top_predictions():
  with open(f'{FLAGS.jobname}/ranking_debug.json', 'r') as json_file:
    scores = json.load(json_file)
  top_pred = scores['order'][:FLAGS.n_top_predictions]
  return top_pred
  
def MF_PAE_matrix():
  all_models_pae = []
  jobname = FLAGS.jobname
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    all_models_pae.append(np.asarray(data['predicted_aligned_error']))
  plot_paes(all_models_pae)
  if FLAGS.action == "save":
    plt.savefig(f"{pred}.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()
    
def MF_plddts():
  all_models_plddt = []
  jobname = FLAGS.jobname
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    all_models_plddt.append(np.asarray(data['plddt']))
  plot_plddts(all_models_plddt)
  if FLAGS.action == "save":
    plt.savefig(f"top_{FLAGS.n_top_predictions}_plddt.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def MF_DM_dual_plddt_PAE(prediction):
  jobname = FLAGS.jobname
  with open(f'{jobname}/result_{prediction}.pkl', "rb") as results_file:
    results = pickle.load(results_file)

  pae_outputs = {}
  pae_outputs["test"] = (results["predicted_aligned_error"], results['max_predicted_aligned_error'])
  pae, max_pae = list(pae_outputs.values())[0]
  plt.figure(figsize=[8 * 2, 6])
  
  plt.subplot(1, 2, 1)
  plt.plot(results['plddt'])
  plt.title('Predicted LDDT')
  plt.xlabel('Residue')
  plt.ylabel('pLDDT')

  plt.subplot(1, 2, 2)
  plt.imshow(pae, vmin=0., vmax=max_pae, cmap='Greens_r')
  plt.colorbar(fraction=0.046, pad=0.04)
  plt.title('Predicted Aligned Error')
  plt.xlabel('Scored residue')
  plt.ylabel('Aligned residue')
  
  if FLAGS.action == "save":
    plt.savefig(f"{prediction}_plddt_PAE.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def call_dual():
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    MF_DM_dual_plddt_PAE(pred)
  
def MF_indiv_plddt():
  jobname = FLAGS.jobname
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    plot_confidence(data['plddt'])
    if FLAGS.action == "save":
      plt.savefig(f"{pred}_plddt_PAE.png")
      plt.close()
    if FLAGS.action == "show":
      plt.show()
            
def MF_coverage():
  jobname = FLAGS.jobname
  with open(f'{jobname}/features.pkl', 'rb') as f:
    data = pickle.load(f)
  plot_msa_v2(data)
  plt.show()
  

def main(argv):
  if "coverage" in FLAGS.chosen_plots:
    MF_coverage()
  if FLAGS.plot_type == "one_for_all":
    MF_PAE_matrix()
    MF_plddts()
  elif FLAGS.plot_type == "for_each":
    plot_for_each()
      



if __name__ == "__main__":
  """ 
  These functions are based on the following scripts from ColabFold repository:
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py
  """
  #deepmind_plddt()
  app.run(main)
  
