#!/usr/bin/env python3
"""MassiveFold workspace installer."""

import glob
import json
import os
import shutil
import socket

from massivefold import engines

TOOLS = ("AFmassive", "AlphaFold3", "ColabFold")

def package_root():
  return os.path.dirname(os.path.abspath(__file__))

def resolve_path(path_value):
  if "$" in path_value:
    return path_value
  return os.path.abspath(os.path.expanduser(path_value))

def ordered_massivefold_section(params):
  key_order = [
    "run_massivefold",
    "engine_prefix",
    "data_dir",
    "uniref_database",
    "jobfile_headers_dir",
    "output_dir",
    "logs_dir",
    "input_dir",
    "models_to_use",
    "pkl_format",
  ]
  massivefold = params.get("massivefold", {})
  ordered = {}
  for key in key_order:
    if key in massivefold:
      ordered[key] = massivefold[key]
  for key in massivefold:
    if key not in ordered:
      ordered[key] = massivefold[key]
  params["massivefold"] = ordered
  return params

def patch_cluster_params(params, data_dir=None, engine_prefix=None):
  massivefold = params.setdefault("massivefold", {})
  if data_dir is not None:
    massivefold["data_dir"] = resolve_path(data_dir)
  if engine_prefix is not None:
    massivefold["engine_prefix"] = engine_prefix
  # the templates are found in the installed package
  massivefold.pop("jobfile_templates_dir", None)
  massivefold.pop("scripts_dir", None)
  return ordered_massivefold_section(params)

def patch_workspace_params(params, root_dir):
  massivefold = params.setdefault("massivefold", {})
  massivefold["jobfile_headers_dir"] = os.path.abspath(os.path.join(root_dir, "headers"))
  massivefold["output_dir"] = os.path.abspath(os.path.join(root_dir, "output"))
  massivefold["logs_dir"] = os.path.abspath(os.path.join(root_dir, "log"))
  massivefold["input_dir"] = os.path.abspath(os.path.join(root_dir, "input"))
  return ordered_massivefold_section(params)

def read_json(path):
  with open(path, "r", encoding="utf-8") as handle:
    return json.load(handle)

def write_json(path, payload):
  with open(path, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=4)

def param_name(tool, host_is_jeanzay=False):
  return f"jeanzay_{tool}_params.json" if host_is_jeanzay else f"{tool}_params.json"

def site_defaults_dir():
  """Database paths and engine launchers kept from one install to the next (never in the package)."""
  config_home = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
  return os.path.join(config_home, "massivefold", "site_defaults")

def resolve_site_default_path(tool, host_is_jeanzay=False):
  """Site default of the tool if a previous install wrote one, else the packaged parameters."""
  site_path = os.path.join(site_defaults_dir(), param_name(tool, host_is_jeanzay))
  if os.path.isfile(site_path):
    return site_path
  return os.path.join(package_root(), "parallelization", param_name(tool, host_is_jeanzay))

def parse_tool_from_param_basename(basename):
  host_is_jeanzay = basename.startswith("jeanzay_")
  if host_is_jeanzay:
    basename = basename[len("jeanzay_"):]
  suffix = "_params.json"
  if not basename.endswith(suffix):
    raise ValueError(f"Unexpected parameter file name: {basename}")
  return basename[: -len(suffix)], host_is_jeanzay

def write_site_default(tool, params, host_is_jeanzay=False):
  os.makedirs(site_defaults_dir(), exist_ok=True)
  site_path = os.path.join(site_defaults_dir(), param_name(tool, host_is_jeanzay))
  write_json(site_path, params)
  return site_path

def install_root(path_value):
  return os.path.abspath(os.path.expanduser(path_value))

def create_tree(root_dir):
  os.makedirs(os.path.join(root_dir, "input"), exist_ok=True)
  os.makedirs(os.path.join(root_dir, "output"), exist_ok=True)
  os.makedirs(os.path.join(root_dir, "log"), exist_ok=True)

def copy_examples(root_dir):
  filenames = ("H1140.fasta", "multirun_setup.json")
  directories = ("input", "./")
  for name, directory in zip(filenames, directories):
    source = os.path.join(package_root(), "examples", name)
    destination = os.path.join(root_dir, directory, name)
    if os.path.exists(source):
      shutil.copy2(source, destination)
    else:
      print(f"Warning: missing packaged example '{source}'.")

def copy_headers(root_dir, host_is_jeanzay):
  """alignment, jobarray and post_treatment headers: Jean Zay ones on Jean Zay, generic ones elsewhere.
  Existing headers are kept, they may have been adapted to the cluster."""
  source = os.path.join(package_root(), "parallelization", "headers")
  destination = os.path.join(root_dir, "headers")
  os.makedirs(destination, exist_ok=True)
  for job_type in ("alignment", "jobarray", "post_treatment"):
    header = f"example_header_{job_type}_jeanzay.slurm" if host_is_jeanzay else f"{job_type}.slurm"
    target = os.path.join(destination, f"{job_type}.slurm")
    if not os.path.exists(target):
      shutil.copy2(os.path.join(source, header), target)
  for example in glob.glob(os.path.join(source, "example_header_*.slurm")):
    shutil.copy2(example, os.path.join(destination, os.path.basename(example)))

def configure_site_defaults(host_is_jeanzay, data_dirs, launchers):
  """Record the given database paths and launchers; the tools not given keep their previous ones."""
  for tool in TOOLS:
    params = read_json(resolve_site_default_path(tool, host_is_jeanzay))
    params = patch_cluster_params(params, data_dir=data_dirs.get(tool) or None, engine_prefix=launchers.get(tool))
    write_site_default(tool, params, host_is_jeanzay)

def build_param_file(tool, destination, root_dir, host_is_jeanzay=False):
  params = read_json(resolve_site_default_path(tool, host_is_jeanzay))
  params = patch_cluster_params(params)
  params = patch_workspace_params(params, root_dir)
  write_json(destination, params)

def link_jeanzay_af3_databases():
  af3_datadir = os.path.join(os.path.expanduser("~"), "af3_datadir")
  os.makedirs(af3_datadir, exist_ok=True)
  for database in glob.glob(os.path.join(os.environ.get("DSDIR", ""), "Alphafold3", "*")):
    link = os.path.join(af3_datadir, os.path.basename(database))
    if not os.path.lexists(link):
      os.symlink(database, link)

def validate_db_paths(args, host_is_jeanzay):
  if host_is_jeanzay:
    return None

  provided = False
  for path in (args.alphafold_databases, args.alphafold3_databases, args.colabfold_databases):
    if path:
      provided = True
      if not os.path.isdir(path):
        return f"{path} doesn't exists"
  if not provided and not args.only_envs:
    return "At least one of --alphafold-db, --alphafold3-db or --colabfold-db is required."
  return None

def is_jeanzay_host():
  return socket.gethostname()[:8] == "jean-zay"

def install_workspace(args):
  host_is_jeanzay = is_jeanzay_host()
  error = validate_db_paths(args, host_is_jeanzay)
  if error:
    print(error)
    return 1

  data_dirs = {
    "AFmassive": args.alphafold_databases,
    "AlphaFold3": args.alphafold3_databases,
    "ColabFold": args.colabfold_databases,
  }
  # the engines of the given databases, all of them on Jean Zay or with only --only-envs
  tools = [ tool for tool in TOOLS if data_dirs[tool] ]
  if host_is_jeanzay or not tools:
    tools = list(TOOLS)

  root_dir = install_root(args.install_path)
  launchers = {}
  if args.no_env:
    print("No engine environment installed: the engines are run from the PATH (engine_prefix is left as is).")
  else:
    engines_dir = args.engines_dir or os.path.join(root_dir, "engines")
    launchers = engines.install_engines(engines_dir, tools, reuse_conda=not args.ignore_conda_envs)

  configure_site_defaults(host_is_jeanzay, data_dirs, launchers)

  if args.only_envs:
    print(f"Engine environments installed for {', '.join(tools)}; workspace file creation skipped (--only-envs).")
    return 0

  create_tree(root_dir)
  copy_examples(root_dir)
  copy_headers(root_dir, host_is_jeanzay)

  if host_is_jeanzay:
    print("Currently on Jean Zay cluster, using prebuilt headers and json parameter files.")
    link_jeanzay_af3_databases()
  for tool in tools:
    build_param_file(tool, os.path.join(root_dir, f"{tool}_params.json"), root_dir, host_is_jeanzay)
  return 0
