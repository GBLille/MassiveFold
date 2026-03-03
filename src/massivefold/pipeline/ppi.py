"""Implementation of `massivefold screen` pipeline."""

from argparse import Namespace

from .screen import screening_pipeline
from .run import run_pipeline
from massivefold.parallelization.unifier import ppi_create_input

def run_ppi_pipeline_internal(args, forwarded_args, run_args, screen_args, scheduler):
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

  # update temporary args partially parsed from cli.py module
  args_updates = {}
  # small molecules present is a screen call
  if context_file:
    args_updates["ligands"] = context_file
    for ppi in ppi_sequences:
      args_updates["sequence"] = ppi
      vars(screen_args).update(args_updates)
      screening_pipeline(
        screen_args,
        forwarded_args,
        scheduler
      )
  # no small molecules is a simple run call
  else:
    args_updates["run_name"] = "PPI"
    for ppi in ppi_sequences:
      args_updates["sequence"] = ppi
      vars(run_args).update(args_updates)
      print(run_args)
      run_pipeline(
        run_args,
        forwarded_args,
        scheduler
      )
  return 0

def ppi_pipeline(args, forwarded_args, run_args, screen_args, scheduler):
  try:
    return run_ppi_pipeline_internal(args, forwarded_args, run_args, screen_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
