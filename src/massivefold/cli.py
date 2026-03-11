#!/usr/bin/env python3
"""CLI for MassiveFold."""

import argparse
import sys

from massivefold.pipeline import run_pipeline
from massivefold.pipeline import multirun_pipeline
from massivefold.pipeline import screening_pipeline
from massivefold.pipeline import ppi_pipeline
from massivefold.install import install_workspace
from massivefold.scheduling import resolve_scheduler

def add_run_arguments(run_parser):
  # required arguments
  run_required = run_parser.add_argument_group("Required arguments")
  run_required.add_argument("-s", "--sequence", dest="sequence", required=True,
                              help="Path of the sequence(s) to infer, should be a 'fasta' file.")
  run_required.add_argument("-r", "--run", dest="run_name", required=True,
                              help="Name chosen for the run to organize in outputs.")
  run_required.add_argument("-f", "--parameters", dest="parameters", required=True,
                              help="Json file's path containing the parameters used for this run.")

  # optional arguments
  run_optional = run_parser.add_argument_group("Optional arguments")
  run_optional.add_argument("-p", "--predictions_per_model", dest="predictions_per_model", type=int, default=5,
                              help="Number of predictions (default: %(default)s) computed for each neural network model."
                              " If used with -t AlphaFold3, -p is the number of seeds used. Each seed will have m"
                              " samples predicted. The number of sample set m is set in the AlphaFold3_params.json file."
                              " In total, with -p n, you will have m*n predictions computed.")
  run_optional.add_argument("-b", "--batch_size", dest="batch_size", type=int, default=25,
                              help="Number of predictions per batch (default: %(default)s). For AlphaFold3, it corresponds to"
                              " number of seeds should not be higher than -p.")
  run_optional.add_argument("-j", "--jobid", dest="jobid", 
                            help="Jobid of an alignment job to wait for inference, skips the alignments.")
  run_optional.add_argument("-o", "--only_msas", dest="only_msas", action="store_true",
                              help="Only compute alignments, the first step of MassiveFold."
                              " Overwrite MSAs directory by forcing re-computation.")
  run_optional.add_argument("-c", "--calibrate", dest="calibrate", action="store_true",
                              help="Calibrate --batch_size value. Searches from the previous runs for the same 'fasta'"
                              " path given in --sequence and uses the longest prediction time found to compute the maximal"
                              " number of predictions per batch. This maximal number depends on the total time given by"
                              " --wall_time.")
  run_optional.add_argument("-C", "--calibration_from", dest="calibration_from",
                              help="Path of a previous run to calibrate the batch size from (see --calibrate).")
  run_optional.add_argument("-w", "--wall_time", dest="wall_time", type=float, default=20,
                              help="Total time in hour (default: %(default)s) available for calibration computations.")
  run_optional.add_argument("-m", "--msas_precomputed", dest="msas_precomputed", 
                              help="Path to directory that contains computed msas.")
  run_optional.add_argument("-n", "--top_n_model", dest="top_n_model",
                              help="Uses the n neural network models with best ranking confidence from this run's path.")
  run_optional.add_argument("-a", "--recompute_msas", dest="recompute_msas", action="store_true",
                              help="Purges previous alignment step and recomputes msas.")
  run_optional.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: %(default)s)",
  )

def add_multirun_arguments(multirun_parser):
  # required arguments
  multirun_required = multirun_parser.add_argument_group("Required arguments")
  multirun_required.add_argument("-s", "--sequence", dest="sequence", required=True,
                                    help="Path of the fasta file containing sequence(s) used screening.")
  multirun_required.add_argument("--setup", dest="setup", required=True,
                                    help="Csv file containing the list of runs with their parameters.")
  multirun_required.add_argument("-f", "--parameters", dest="parameters", required=True,
                              help="Json file's path containing the parameters used for this run.")

  # optional arguments
  multirun_optional = multirun_parser.add_argument_group("Optional arguments")
  multirun_optional.add_argument("-p", "--predictions_per_model", dest="predictions_per_model", type=int, default=1,
                                    help="(default: %(default)s) Number of seed used with AlphaFold3."
                                    " Each seed will have 5 samples predicted. In total, with -p n,"
                                    " you will have 5n predictions computed (5 predictions with default params).")
  multirun_optional.add_argument("-m", "--msas_precomputed", dest="msas_precomputed",
                                    help="Path to directory that contains computed msas.")
  multirun_optional.add_argument("-o", "--only_msas", dest="only_msas", action="store_true",
                                    help="Only compute alignments, the first step of MassiveFold.")
  multirun_optional.add_argument("-j", "--jobid", dest="jobid",
                                    help="Jobid of an alignment job to wait for inference, skips the alignments.")
  multirun_optional.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: %(default)s)",
  )

def add_install_arguments(install_parser):
  install_parser.add_argument("--alphafold-db", dest="alphafold_databases", default="", help="Path to AlphaFold2 database.")
  install_parser.add_argument("--alphafold3-db", dest="alphafold3_databases", default="", help="Path to AlphaFold3 database.")
  install_parser.add_argument("--colabfold-db", dest="colabfold_databases", default="", help="Path to ColabFold database.")
  install_parser.add_argument("--install-path", dest="install_path", default="massivefold_runs",
                                help="Where to install MassiveFold files.")
  install_parser.add_argument("--no-env", dest="no_env", action="store_true",
                                help="No environments installation but only files and parameters. At least one of --alphafold-db"
                                " or --colabfold-db is required with this option.")
  install_parser.add_argument("--only-envs", dest="only_envs", action="store_true",
                                help="Only install the environments (other arguments are not used.")

def add_screening_arguments(screening_parser):
  # required arguments
  screening_required = screening_parser.add_argument_group("Required arguments")
  screening_required.add_argument("-s", "--sequence", dest="sequence", required=True,
                                    help="Path of the fasta file containing sequence(s) used for screening.")
  screening_required.add_argument("-l", "--ligands", dest="ligands", required=True,
                                    help="Csv file containing the list of ligands to use for screening.")
  screening_required.add_argument("-f", "--parameters", dest="parameters", required=True,
                                    help="Json file's path containing the parameters used for the screening.")

  # optional arguments
  screening_optional = screening_parser.add_argument_group("Optional arguments")
  screening_optional.add_argument("-p", "--predictions_per_model", dest="predictions_per_model", type=int, default=1,
                                    help="(default: %(default)s) Number of seed used with AlphaFold3."
                                    " Each seed will have 5 samples predicted. In total, with -p n,"
                                    " you will have 5n predictions computed (5 predictions with default params).")
  screening_optional.add_argument("-m", "--msas_precomputed", dest="msas_precomputed",
                                    help="Path to directory that contains computed msas.")
  screening_optional.add_argument("-o", "--only_msas", dest="only_msas", action="store_true",
                                    help="Only compute alignments, the first step of MassiveFold.")
  screening_optional.add_argument("-j", "--jobid", dest="jobid",
                                    help="Jobid of an alignment job to wait for inference, skips the alignments.")
  screening_optional.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: %(default)s)",
  )

def add_ppi_arguments(ppi_parser):
  # required arguments
  ppi_required = ppi_parser.add_argument_group('Required arguments')
  ppi_required.add_argument("--receptors", dest="receptors", required=True,
                              help="Path to CSV file containing fasta file paths to proteic sequence(s) used as receptors.")
  ppi_required.add_argument("--ligands", dest="ligands", required=True,
                              help="Path to CSV file containing fasta file paths to proteic sequence(s) used as ligands.")
  ppi_required.add_argument("-f", "--parameters", dest="parameters", required=True,
                              help="Json file's path containing the parameters used for the screening.")
  # optional arguments
  ppi_optional = ppi_parser.add_argument_group('Optional arguments')
  ppi_optional.add_argument("--context", dest="context",
                              help="Path to CSV file containing molecules that are used as context (substrate, ions, ...) in the PPI simulations.")
  ppi_optional.add_argument("-p", "--predictions_per_model", dest="predictions_per_model", type=int, default=1,
                              help="Number of predictions to generate. For AlphaFold3, this sets the number of seeds."
                              " For AlphaFold2, it sets the number of predictions per model (default: %(default)s)")
  ppi_optional.add_argument("-m", "--msas_precomputed", dest="msas_precomputed",
                              help="Path to directory that contains computed msas.")
  ppi_optional.add_argument("-o", "--only_msas", dest="only_msas", action="store_true",
                              help="Only compute alignments, the first step of MassiveFold.")
  ppi_optional.add_argument("-j", "--jobid", dest="jobid",
                              help="Jobid of an alignment job to wait for inference, skips the alignments.")
  ppi_optional.add_argument(
    "--scheduler",
    dest="scheduler",
    default="auto",
    choices=["auto", "slurm", "local"],
    help="Scheduler selector (default: auto) (default: %(default)s)",
  )

def build_parser():
  parser = argparse.ArgumentParser(prog="massivefold")
  subparsers = parser.add_subparsers(dest="command")

  run_parser = subparsers.add_parser("run", help="Run MassiveFold")
  add_run_arguments(run_parser)

  multirun_parser = subparsers.add_parser("multirun", help="Run MassiveFold")
  add_multirun_arguments(multirun_parser)

  screening_parser = subparsers.add_parser("screen", help="Run MassiveFold screening")
  add_screening_arguments(screening_parser)

  ppi_parser = subparsers.add_parser("ppi", help="Run MassiveFold PPI screening")
  add_ppi_arguments(ppi_parser)

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

def dispatch_multirun(args, forwarded_args, scheduler):
  print(f"Selected scheduler: {scheduler['name']}")
  try:
    return multirun_pipeline(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    return 1

def dispatch_screen(args, forwarded_args, scheduler):
  print(f"Selected scheduler: {scheduler['name']}")
  try:
    return screening_pipeline(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    return 1

def dispatch_ppi(args, forwarded_args, scheduler):
  print(f"Selected scheduler: {scheduler['name']}")
  try:
    return ppi_pipeline(args, forwarded_args, scheduler)
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

  if args.command in ["run", "multirun", "screen", "ppi"]:
    scheduler, status = resolve_selected_scheduler(args)
    if status != 0:
      return status

    if args.command == "run":
      return dispatch_run(args, unknown, scheduler)
    if args.command == "multirun":
      return dispatch_multirun(args, unknown, scheduler)
    elif args.command == "screen":
      return dispatch_screen(args, unknown, scheduler)
    elif args.command == "ppi":
      return dispatch_ppi(args, unknown, scheduler)

  if args.command == "install":
    if unknown:
      print("Ignoring extra arguments:", " ".join(unknown))
    return dispatch_install(args)

  parser.print_help()
  return 1

if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
