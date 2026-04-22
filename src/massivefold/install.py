#!/usr/bin/env python3
"""MassiveFold workspace installer."""

import json
import os
import shutil
import socket

def package_root():
  return os.path.dirname(os.path.abspath(__file__))

def site_defaults_root():
  return os.path.join(package_root(), "parallelization", "site_defaults")

def tool_param_name(tool, host_is_jeanzay=False):
  base_name = f"{tool}_params.json"
  if host_is_jeanzay:
    return f"jeanzay_{base_name}"
  return base_name

def site_default_path(tool):
  return os.path.join(site_defaults_root(), tool_param_name(tool))

def resolve_path(path_value):
  if "$" in path_value:
    return path_value
  return os.path.abspath(os.path.expanduser(path_value))

def ordered_massivefold_section(params):
  key_order = [
    "run_massivefold",
    "data_dir",
    "uniref_database",
    "jobfile_headers_dir",
    "jobfile_templates_dir",
    "scripts_dir",
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

def patch_cluster_params(params, data_dir=None):
  pkg_root = package_root()
  massivefold = params.setdefault("massivefold", {})

  if data_dir is not None:
    massivefold["data_dir"] = resolve_path(data_dir)
  massivefold["jobfile_templates_dir"] = os.path.abspath(os.path.join(pkg_root, "parallelization", "templates"))
  massivefold["scripts_dir"] = os.path.abspath(os.path.join(pkg_root, "parallelization"))
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

def install_root(path_value):
  return os.path.abspath(os.path.expanduser(path_value))

def create_tree(root_dir):
  os.makedirs(os.path.join(root_dir, "input"), exist_ok=True)
  os.makedirs(os.path.join(root_dir, "output"), exist_ok=True)
  os.makedirs(os.path.join(root_dir, "log"), exist_ok=True)

def copy_example_fasta(root_dir):
  source = os.path.join(package_root(), "examples", "H1140.fasta")
  destination = os.path.join(root_dir, "input", "H1140.fasta")
  if os.path.exists(source):
    shutil.copy2(source, destination)
  else:
    print(f"Warning: missing packaged example '{source}'.")

def copy_headers(root_dir):
  source = os.path.join(package_root(), "parallelization", "headers")
  destination = os.path.join(root_dir, "headers")
  os.makedirs(destination, exist_ok=True)
  for name in os.listdir(source):
    if name.endswith(".slurm"):
      shutil.copy2(os.path.join(source, name), os.path.join(destination, name))

def configure_site_defaults(host_is_jeanzay, alphafold_db=None, alphafold3_db=None, colabfold_db=None):
  os.makedirs(site_defaults_root(), exist_ok=True)
  data_dirs = {
    "AFmassive": alphafold_db,
    "AlphaFold3": alphafold3_db,
    "ColabFold": colabfold_db,
  }
  for tool, data_dir in data_dirs.items():
    params = read_json(os.path.join(
      package_root(),
      "parallelization",
      tool_param_name(tool, host_is_jeanzay)
    ))
    params = patch_cluster_params(params, data_dir=data_dir)
    write_json(site_default_path(tool), params)

def preferred_param_source(tool, host_is_jeanzay=False):
  candidate = site_default_path(tool)
  if os.path.exists(candidate):
    return candidate
  return os.path.join(package_root(), "parallelization", tool_param_name(tool, host_is_jeanzay))

def build_param_file(tool, destination, root_dir, data_dir=None, host_is_jeanzay=False):
  params = read_json(preferred_param_source(tool, host_is_jeanzay))
  if not os.path.exists(site_default_path(tool)):
    params = patch_cluster_params(params, data_dir=data_dir)
  params = patch_workspace_params(params, root_dir)
  write_json(destination, params)

def install_jeanzay(root_dir):
  build_param_file(
    "AFmassive",
    os.path.join(root_dir, "AFmassive_params.json"),
    root_dir,
    host_is_jeanzay=True,
  )
  build_param_file(
    "AlphaFold3",
    os.path.join(root_dir, "AlphaFold3_params.json"),
    root_dir,
    host_is_jeanzay=True,
  )
  build_param_file(
    "ColabFold",
    os.path.join(root_dir, "ColabFold_params.json"),
    root_dir,
    host_is_jeanzay=True,
  )

  headers = os.path.join(root_dir, "headers")
  source_to_target = {
    "example_header_alignment_jeanzay.slurm": "alignment.slurm",
    "example_header_jobarray_jeanzay.slurm": "jobarray.slurm",
    "example_header_post_treatment_jeanzay.slurm": "post_treatment.slurm",
  }
  for source_name, target_name in source_to_target.items():
    source = os.path.join(headers, source_name)
    target = os.path.join(headers, target_name)
    if os.path.exists(source) and not os.path.exists(target):
      os.rename(source, target)

def install_standard(root_dir, alphafold_db=None, alphafold3_db=None, colabfold_db=None):
  if alphafold_db:
    build_param_file(
      "AFmassive",
      os.path.join(root_dir, "AFmassive_params.json"),
      root_dir,
      data_dir=alphafold_db,
    )
  if alphafold3_db:
    build_param_file(
      "AlphaFold3",
      os.path.join(root_dir, "AlphaFold3_params.json"),
      root_dir,
      data_dir=alphafold3_db,
    )
  if colabfold_db:
    build_param_file(
      "ColabFold",
      os.path.join(root_dir, "ColabFold_params.json"),
      root_dir,
      data_dir=colabfold_db,
    )

def validate_db_paths(args, host_is_jeanzay):
  if host_is_jeanzay:
    return None

  provided = False
  if args.alphafold_databases:
    provided = True
    if not os.path.isdir(args.alphafold_databases):
      return f"{args.alphafold_databases} doesn't exists"
  if args.alphafold3_databases:
    provided = True
    if not os.path.isdir(args.alphafold3_databases):
      return f"{args.alphafold3_databases} doesn't exists"
  if args.colabfold_databases:
    provided = True
    if not os.path.isdir(args.colabfold_databases):
      return f"{args.colabfold_databases} doesn't exists"
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

  configure_site_defaults(
    host_is_jeanzay,
    alphafold_db=args.alphafold_databases,
    alphafold3_db=args.alphafold3_databases,
    colabfold_db=args.colabfold_databases,
  )

  if args.only_envs:
    print("Configured internal MassiveFold defaults and skipped workspace file creation because --only-envs is set.")
    return 0

  root_dir = install_root(args.install_path)
  create_tree(root_dir)
  copy_example_fasta(root_dir)
  copy_headers(root_dir)

  if host_is_jeanzay:
    print("Currently on Jean Zay cluster, using prebuilt headers and json parameter files.")
    install_jeanzay(root_dir)
    return 0

  install_standard(
    root_dir,
    alphafold_db=args.alphafold_databases,
    alphafold3_db=args.alphafold3_databases,
    colabfold_db=args.colabfold_databases,
  )
  return 0
