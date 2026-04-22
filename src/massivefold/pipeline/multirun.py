"""Implementation of `massivefold multirun` pipeline."""

from argparse import Namespace
import copy
import json
import os

from .run import run_pipeline
from massivefold.parallelization.unifier import get_multirun_runs

def run_item_args(args, run_name, run_args, parameters_file):
  allowed_args = {
    "predictions_per_model",
    "batch_size",
    "jobid",
    "only_msas",
    "calibrate",
    "calibration_from",
    "wall_time",
    "msas_precomputed",
    "top_n_model",
    "recompute_msas",
  }
  unknown_args = [ key for key in run_args if key not in allowed_args ]
  if unknown_args:
    raise ValueError(
      f"Unknown run arg(s) for '{run_name}': {', '.join(unknown_args)}"
    )
  return Namespace(
    sequence=args.sequence,
    run_name=run_name,
    parameters=parameters_file,
    predictions_per_model=run_args.get("predictions_per_model", 1),
    batch_size=run_args.get("batch_size", 25),
    jobid=run_args.get("jobid", None),
    only_msas=run_args.get("only_msas", False),
    calibrate=run_args.get("calibrate", False),
    calibration_from=run_args.get("calibration_from", None),
    wall_time=run_args.get("wall_time", 20),
    msas_precomputed=run_args.get("msas_precomputed", None),
    top_n_model=run_args.get("top_n_model", None),
    recompute_msas=run_args.get("recompute_msas", False),
    scheduler=args.scheduler,
  )

def multirun_pipeline_internal(args, forwarded_args, scheduler):
  setup_json = args.setup
  runs = get_multirun_runs(setup_json)

  status = 0
  setup_dir = os.path.dirname(os.path.abspath(setup_json))
  for run in runs:
    sequence_parameters = copy.deepcopy(run["parameters"])
    file_id = os.urandom(6).hex()
    internal_parameters_file = os.path.join(
      setup_dir,
      f"{file_id}.json"
    )
    try:
      json.dump(sequence_parameters, open(internal_parameters_file, 'w'), indent=4)
      pipeline_args = run_item_args(
        args,
        run["name"],
        run["args"],
        internal_parameters_file
      )
      run_status = run_pipeline(
        pipeline_args,
        forwarded_args,
        scheduler
      )
      status = max(status, run_status)
    finally:
      if os.path.exists(internal_parameters_file):
        os.remove(internal_parameters_file)

  return status

def multirun_pipeline(args, forwarded_args, scheduler):
  try:
    return multirun_pipeline_internal(args, forwarded_args, scheduler)
  except (RuntimeError, ValueError) as error:
    print(error)
    print("Exiting.")
    return 1
