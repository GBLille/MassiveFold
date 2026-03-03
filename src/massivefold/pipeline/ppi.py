"""Implementation of `massivefold screen` pipeline."""

from argparse import Namespace

from .screen import screening_pipeline
from .run import run_pipeline
from massivefold.parallelization.unifier import ppi_create_input

def run_item_args(args, sequence):
  return Namespace(
    sequence=sequence,
    run_name="PPI",
    parameters=args.parameters,
    predictions_per_model=args.predictions_per_model,
    batch_size=25,
    jobid=args.jobid,
    only_msas=args.only_msas,
    calibrate=False,
    calibration_from=None,
    wall_time=20,
    msas_precomputed=args.msas_precomputed,
    top_n_model=None,
    recompute_msas=False,
    scheduler=args.scheduler,
  )

def screening_item_args(args, sequence):
  return Namespace(
    sequence=sequence,
    ligands=args.context,
    parameters=args.parameters,
    predictions_per_model=args.predictions_per_model,
    msas_precomputed=args.msas_precomputed,
    only_msas=args.only_msas,
    jobid=args.jobid,
    scheduler=args.scheduler,
  )

def run_ppi_pipeline_internal(args, forwarded_args, scheduler):
  receptors_file = args.receptors
  ligands_file = args.ligands
  context_file = args.context
  parameters_file = args.parameters

  ppi_inputs = ppi_create_input(
    receptors_file,
    ligands_file,
    context_file,
    parameters_file
  )

  ppi_sequences = ppi_inputs["ppi"].tolist()

  # small molecules present is a screen call
  if context_file:
    for ppi in ppi_sequences:
      screen_args = screening_item_args(args, ppi)
      screening_pipeline(
        screen_args,
        forwarded_args,
        scheduler
      )
  # no small molecules is a simple run call
  else:
    for ppi in ppi_sequences:
      pipeline_args = run_item_args(args, ppi)
      run_pipeline(
        pipeline_args,
        forwarded_args,
        scheduler
      )
  return 0

def ppi_pipeline(args, forwarded_args, scheduler):
  try:
    return run_ppi_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
