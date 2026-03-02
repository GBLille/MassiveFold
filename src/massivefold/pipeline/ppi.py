"""Implementation of `massivefold screen` pipeline."""

import os
import shutil

from .screen import run_screening_pipeline_internal
from massivefold.parallelization.unifier import ppi_create_input

def run_ppi_pipeline_internal(args, forwared_args, scheduler):
  receptors_file = args.receptors
  ligands_file = args.ligands
  context_file = args.context
  parameters_file = args.parameters

  predictions_per_model = args.predictions_per_model
  msas_precomputed = args.msas_precomputed
  only_msas = args.only_msas
  wait_for_jobid = args.jobid

  a = ppi_create_input(
    receptors_file,
    ligands_file,
    context_file,
    parameters_file
  )

  print(a)

def ppi_pipeline(args, forwarded_args, scheduler):
  try:
    return run_ppi_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
