"""Implementation of `massivefold screen` pipeline."""

from argparse import Namespace
import copy
import json
import os

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
    parameters_file
  )

  ppi_sequences = ppi_inputs["ppi"].tolist()
  ppi_fasta_chains = ppi_inputs["fasta_chains"].tolist()
  parameters = json.load(open(parameters_file, 'r'))
  is_af3 = 'AF3_run' in parameters

  for ppi, fasta_chains in zip(ppi_sequences, ppi_fasta_chains):
    sequence_parameters = copy.deepcopy(parameters)
    # automatically detect and register the fasta chain types in af3 param file
    if is_af3:
      sequence_parameters["AF3_run"]["fasta_chains"] = fasta_chains
    # created combined fasta stored inside input directory
    combined_fasta_dir = os.path.join(sequence_parameters["massivefold"]["input_dir"], 'combined')
    os.makedirs(combined_fasta_dir, exist_ok=True)
    sequence_parameters["massivefold"]["input_dir"] = combined_fasta_dir
    # copy param file in each sequence's log directory
    file_id = os.urandom(6).hex()
    internal_parameters_file = os.path.join(
      os.path.dirname(parameters_file),
      f"{file_id}.json"
    )
    json.dump(sequence_parameters, open(internal_parameters_file, 'w'))
    args.parameters = internal_parameters_file
    # small molecules present is a screen call
    if context_file:
      screen_args = screening_item_args(args, ppi)
      screening_pipeline(
        screen_args,
        forwarded_args,
        scheduler
      )
    # no small molecules is a simple run call
    else:
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
