"""Implementation of `massivefold screen` pipeline."""

from argparse import Namespace
import copy
import json
import os

from .screen import screening_pipeline
from .run import run_pipeline
from .run import detect_tool_code
from massivefold.parallelization.unifier import get_multirun_runs

def run_item_args(args, run):
  return Namespace(
    sequence=args.sequence,
    run_name=run,
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

def multirun_pipeline_internal(args, forwarded_args, scheduler):
  setup_csv = args.setup
  parameters_file = args.parameters
  massivefold_params = json.load(open(parameters_file, 'r'))["massivefold"]
  all_tool, all_runs, all_parameter_files = get_multirun_runs(
    setup_csv, massivefold_params
  )
  for tool, run_name, parameters in zip(all_tool, all_runs, all_parameter_files):
    sequence_parameters = copy.deepcopy(parameters)
    #######
    # copy param file in each sequence's log directory
    file_id = os.urandom(6).hex()
    internal_parameters_file = os.path.join(
      os.path.dirname(parameters_file),
      f"{file_id}.json"
    )
    json.dump(sequence_parameters, open(internal_parameters_file, 'w'))
    args.parameters = internal_parameters_file
    # small molecules present is a screen call
    pipeline_args = run_item_args(args, run_name)
    run_pipeline(
      pipeline_args,
      forwarded_args,
      scheduler
    )
    os.remove(internal_parameters_file)
  return 0

def multirun_pipeline(args, forwarded_args, scheduler):
  try:
    return multirun_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
