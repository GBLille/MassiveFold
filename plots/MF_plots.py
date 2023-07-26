import numpy as np
import os
import sys
import pickle
import json
from absl import flags
from absl import app
import matplotlib.pyplot as plt
from colabfold_plots import plot_msa_v2, plot_plddts, plot_confidence, plot_paes, plot_plddt_legend

FLAGS = flags.FLAGS

flags.DEFINE_string('input_path', None, 'Path to directory were the alphafold output to be plotted lies.')
flags.DEFINE_enum('plot_type', "", ['one_for_all', 'for_each', ""],
                  'Select either a global plot for all specified predictions or one plot for each of them.')
flags.DEFINE_integer('top_n_predictions', 10, 'Specify the number of predictions taken into account for plotting, it will be the n best predictions.')
flags.DEFINE_list('chosen_plots', [], 'Specify the plots you want to get.')
flags.DEFINE_enum('action', "save", ["save", "show"], "Chose to save the plot or show them.")
flags.DEFINE_string('output_path', None, 'Path to directory that will be store the plot, same as the input dir by default.')

PLOT_TYPES = {
  "for_each": ["DM_plddt_PAE", "CF_plddt"],
  "one_for_all": ["CF_PAEs", "CF_plddts"],
  "specific": ["coverage"]
}

def extract_top_predictions():
  with open(f'{FLAGS.input_path}/ranking_debug.json', 'r') as json_file:
    scores = json.load(json_file)
  top_pred = scores['order'][:FLAGS.top_n_predictions]
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

def MF_DM_dual_plddt_PAE(prediction):
  jobname = FLAGS.input_path
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
    plt.savefig(f"{FLAGS.output_path}/{prediction}_plddt_PAE.png")
    print(f"Saved as {prediction}_plddt_PAE.png")
    plt.close()
  if FLAGS.action == "show":
    plt.show()

def call_dual():
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    MF_DM_dual_plddt_PAE(pred)
  
def MF_indiv_plddt():
  jobname = FLAGS.input_path
  preds_to_plot = extract_top_predictions()
  for pred in preds_to_plot:
    with open(f'{jobname}/result_{pred}.pkl', "rb") as pkl_file:
      data = pickle.load(pkl_file)
    plot_confidence(data['plddt'])
    if FLAGS.action == "save":
      plt.savefig(f"{FLAGS.output_path}/{pred}_plddt.png")
      print(f"Saved as {pred}_plddt.png")
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
  if FLAGS.action == "show":
    plt.show()

PLOT_MAP = {
  "DM_plddt_PAE": call_dual,
  "CF_plddt": MF_indiv_plddt,
  "CF_PAEs": CF_PAEs,
  "CF_plddts": CF_plddts
}

def MF_plot():
  # plot names assignment
  chosen = []
  if FLAGS.chosen_plots:
    chosen = FLAGS.chosen_plots
  # if there is no plot names from the plot_type selected, add them by default
  if FLAGS.plot_type and not set(chosen) & set(PLOT_TYPES[FLAGS.plot_type]):
    chosen.extend(PLOT_TYPES[FLAGS.plot_type])

  # plot all the chosen plots
  for plot in chosen:
    if plot not in PLOT_TYPES['specific']:
      PLOT_MAP[plot]()

def main(argv):
  FLAGS.input_path = os.path.realpath(FLAGS.input_path)
  if not FLAGS.chosen_plots and not FLAGS.plot_type:
    raise ValueError(f"Chose either a plot type to have the default plots associated or\
      directly chose the plots you want")
  if not FLAGS.output_path:
    FLAGS.output_path = FLAGS.input_path

  # Plot depending on the user specifications
  if "coverage" in FLAGS.chosen_plots:
    MF_coverage()
  MF_plot()




if __name__ == "__main__":
  """ 
  These functions are based on the following scripts from ColabFold repository and DeepMind colab notebook:
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/plot.py
  https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py
  https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb
  """
  app.run(main)
  
