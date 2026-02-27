#!/usr/bin/env python3
"""CLI for MassiveFold."""

import argparse
import sys

from massivefold.pipeline import run_pipeline
from massivefold.install import install_workspace
from massivefold.scheduling import resolve_scheduler

def add_run_arguments(run_parser):
  run_parser.add_argument("-s", "--sequence", dest="sequence", required=True, help="Path of the input FASTA file")
  run_parser.add_argument("-r", "--run", dest="run_name", required=True, help="Run name")
  run_parser.add_argument("-f", "--parameters", dest="parameters", required=True, help="Path to parameter JSON file")
  run_parser.add_argument("-p", "--predictions_per_model", dest="predictions_per_model", type=int, default=5)
  run_parser.add_argument("-b", "--batch_size", dest="batch_size", type=int, default=25)
  run_parser.add_argument("-c", "--calibrate", dest="calibrate", action="store_true")
  run_parser.add_argument("-C", "--calibration_from", dest="calibration_from")
  run_parser.add_argument("-w", "--wall_time", dest="wall_time", type=float, default=20)
  run_parser.add_argument("-m", "--msas_precomputed", dest="msas_precomputed")
  run_parser.add_argument("-n", "--top_n_model", dest="top_n_model")
  run_parser.add_argument("-a", "--recompute_msas", dest="recompute_msas", action="store_true")
  run_parser.add_argument("-o", "--only_msas", dest="only_msas", action="store_true")
  run_parser.add_argument("-j", "--jobid", dest="jobid")
  run_parser.add_argument("-t", "--tool", dest="tool")
  run_parser.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: auto)",
  )

def add_install_arguments(install_parser):
  install_parser.add_argument("--alphafold-db", dest="alphafold_databases", default="")
  install_parser.add_argument("--alphafold3-db", dest="alphafold3_databases", default="")
  install_parser.add_argument("--colabfold-db", dest="colabfold_databases", default="")
  install_parser.add_argument("--install-path", dest="install_path", default="massivefold_runs")
  install_parser.add_argument("--no-env", dest="no_env", action="store_true")
  install_parser.add_argument("--only-envs", dest="only_envs", action="store_true")

def build_parser():
  parser = argparse.ArgumentParser(prog="massivefold")
  subparsers = parser.add_subparsers(dest="command")

  run_parser = subparsers.add_parser("run", help="Run MassiveFold")
  add_run_arguments(run_parser)

  screening_parser = subparsers.add_parser("screening", help="Run MassiveFold screening")
  screening_parser.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: auto)",
  )

  install_parser = subparsers.add_parser("install", help="Create MassiveFold file architecture")
  add_install_arguments(install_parser)

  return parser

def dispatch_run(args, forwarded_args, scheduler):
  print(f"Selected scheduler: {scheduler['name']}")
  try:
    return run_pipeline(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    return 1

def dispatch_screening(args, forwarded_args, scheduler):
  print(f"Selected scheduler: {scheduler['name']}")
  try:
    return scheduler["screening"](args, forwarded_args)
  except RuntimeError as error:
    print(error)
    return 1

def dispatch_install(args):
  return install_workspace(args)

def resolve_selected_scheduler(args):
  scheduler, error = resolve_scheduler(args.scheduler)
  if error:
    print(error)
    return None, 2
  return scheduler, 0

def main(argv=None):
  parser = build_parser()
  args, unknown = parser.parse_known_args(argv)

  if not args.command:
    parser.print_help()
    return 0

  if args.command in ["run", "screening"]:
    scheduler, status = resolve_selected_scheduler(args)
    if status != 0:
      return status

    if args.command == "run":
      return dispatch_run(args, unknown, scheduler)
    return dispatch_screening(args, unknown, scheduler)

  if args.command == "install":
    if unknown:
      print("Ignoring extra arguments:", " ".join(unknown))
    return dispatch_install(args)

  parser.print_help()
  return 1

if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
