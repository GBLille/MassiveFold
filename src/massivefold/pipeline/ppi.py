"""Implementation of `massivefold screen` pipeline."""

from argparse import Namespace

from .screen import screening_pipeline
from .run import run_pipeline
from massivefold.parallelization.unifier import ppi_create_input

def run_ppi_pipeline_internal(args, forwarded_args, scheduler):
  receptors_file = args.receptors
  ligands_file = args.ligands
  context_file = args.context
  parameters_file = args.parameters

  predictions_per_model = args.predictions_per_model
  msas_precomputed = args.msas_precomputed
  only_msas = args.only_msas
  wait_for_jobid = args.jobid

  ppi_inputs = ppi_create_input(
    receptors_file,
    ligands_file,
    context_file,
    parameters_file
  )

  ppi_sequences = ppi_inputs["ppi"].tolist()
  base_args = {
    "parameters": parameters_file,
    "predictions_per_model": predictions_per_model,
    "msas_precomputed": msas_precomputed,
    "only_msas": only_msas,
    "wait_for_jobid": wait_for_jobid
  }
  for ppi in ppi_sequences:
    args_ppi = base_args.copy()
    args_ppi["sequence"] = ppi
    if context_file:
      args_ppi["ligands"] = context_file
      args = Namespace(**args_ppi)
      screening_pipeline(args, forwarded_args, scheduler)
    else:
      args_ppi["run"] = "PPI"
      args = Namespace(**args_ppi)
      run_pipeline(args, forwarded_args, scheduler)

  return 0

def ppi_pipeline(args, forwarded_args, scheduler):
  try:
    return run_ppi_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
