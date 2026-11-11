"""
Each engine runs from its own pixi environment, (see pixi.toml and pixi.lock).
Installed when `massivefold install` is run in an engines directory with one 
launcher per engine, used as `engine_prefix` in the jobfiles (`<launcher> <command> <args>`).

If the engine's conda environment (`mf-afmassive-<version>`, ...) is already installed,
it will be used instead.
"""
import glob
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import tarfile
import tempfile
import urllib.request

PIXI_VERSION = "0.82.1"
PIXI_SHA256 = {
  "x86_64": "7ac72af79c50653bfcfddb2473b6e9eff2654affe0057d19fe21b0d8556a7fe2",
}

AFMASSIVE_VERSION = "1.1.11"
COLABFOLD_VERSION = "1.6.3"
ALPHAFOLD3_VERSION = "1.2.0"

# scripts to download
ENGINES = {
  "AFmassive": {
    "environment": "afmassive",
    "conda_environment": f"mf-afmassive-{AFMASSIVE_VERSION}",
    "command": "run_AFmassive.py",
    "distribution": None,
    "scripts": {
      "run_AFmassive.py": f"https://raw.githubusercontent.com/GBLille/AFmassive/v{AFMASSIVE_VERSION}/run_AFmassive.py",
    },
    "package_files": {
      "alphafold/common/stereo_chemical_props.txt":
        f"https://raw.githubusercontent.com/GBLille/AFmassive/v{AFMASSIVE_VERSION}/alphafold/common/stereo_chemical_props.txt",
    },
    "post_install": [],
  },
  "ColabFold": {
    "environment": "colabfold",
    "conda_environment": f"mf-colabfold-{COLABFOLD_VERSION}",
    "command": "colabfold_batch",
    "distribution": ("colabfold", COLABFOLD_VERSION),
    "scripts": {},
    "package_files": {},
    "post_install": [],
  },
  "AlphaFold3": {
    "environment": "alphafold3",
    "conda_environment": f"mf-alphafold3-{ALPHAFOLD3_VERSION}",
    "command": "run_alphafold.py",
    "distribution": ("alphafold3", ALPHAFOLD3_VERSION),
    "scripts": {
      "run_alphafold.py": f"https://raw.githubusercontent.com/GBLille/mf-alphafold3/v{ALPHAFOLD3_VERSION}/run_alphafold.py",
    },
    "package_files": {},
    "post_install": [["build_data"]],
  },
}

def manifest_dir():
  return os.path.dirname(os.path.abspath(__file__))

def launcher_path(engines_dir, tool):
  return os.path.join(engines_dir, "bin", f"run_in_{ENGINES[tool]['environment']}")

def environment_prefix(engines_dir, tool):
  return os.path.join(engines_dir, ".pixi", "envs", ENGINES[tool]["environment"])

def download(url, destination):
  with urllib.request.urlopen(url) as response, open(destination, "wb") as output:
    shutil.copyfileobj(response, output)

def pixi_binary(engines_dir):
  if os.environ.get("MASSIVEFOLD_PIXI"):
    return os.environ["MASSIVEFOLD_PIXI"]
  binary = os.path.join(engines_dir, "bin", "pixi")
  if os.path.isfile(binary):
    return binary
  machine = platform.machine()
  if platform.system() != "Linux" or machine not in PIXI_SHA256:
    raise RuntimeError(f"The engines can only be installed on Linux x86_64, not {platform.system()} {machine}.")
  archive_name = f"pixi-{machine}-unknown-linux-musl.tar.gz"
  url = f"https://github.com/prefix-dev/pixi/releases/download/v{PIXI_VERSION}/{archive_name}"
  print(f"Downloading pixi {PIXI_VERSION}", flush=True)
  os.makedirs(os.path.dirname(binary), exist_ok=True)
  with tempfile.TemporaryDirectory() as tmp:
    archive = os.path.join(tmp, archive_name)
    download(url, archive)
    with open(archive, "rb") as handle:
      digest = hashlib.sha256(handle.read()).hexdigest()
    if digest != PIXI_SHA256[machine]:
      raise RuntimeError(f"Checksum mismatch for {url}: {digest}")
    with tarfile.open(archive) as tar:
      member = tar.getmember("pixi")
      member.name = os.path.basename(binary)
      tar.extract(member, os.path.dirname(binary))
  os.chmod(binary, os.stat(binary).st_mode | stat.S_IXUSR)
  return binary

def pixi_env(engines_dir):
  env = dict(os.environ)
  env.setdefault("CONDA_OVERRIDE_CUDA", "12")
  if not env.get("PIXI_CACHE_DIR") and not env.get("RATTLER_CACHE_DIR"):
    env["PIXI_CACHE_DIR"] = os.path.join(engines_dir, ".cache")
  return env

def run(command, env):
  print(" ".join(command), flush=True)
  subprocess.run(command, check=True, env=env)

def write_launcher(engines_dir, tool, activation_command, env=None):
  """Launcher sourcing the activation script printed by activation_command (pixi shell-hook or conda)."""
  environment = ENGINES[tool]["environment"]
  activation = os.path.join(engines_dir, f"activate_{environment}.sh")
  with open(activation, "w") as handle:
    subprocess.run(activation_command, check=True, env=env, stdout=handle)
  launcher = launcher_path(engines_dir, tool)
  with open(launcher, "w") as handle:
    handle.write(
      "#!/bin/bash\n"
      f"# Runs a command in the '{environment}' environment of the MassiveFold engines\n"
      f"source \"{activation}\"\n"
      "exec \"$@\"\n")
  os.chmod(launcher, 0o755)
  return launcher

def add_engine_files(engines_dir, tool):
  prefix = environment_prefix(engines_dir, tool)
  for name, url in ENGINES[tool]["scripts"].items():
    script = os.path.join(prefix, "bin", name)
    download(url, script)
    with open(script, "r") as handle:
      content = handle.read()
    if not content.startswith("#!"):
      content = "#!/usr/bin/env python\n" + content
      with open(script, "w") as handle:
        handle.write(content)
    os.chmod(script, 0o755)
  for relative_path, url in ENGINES[tool]["package_files"].items():
    site_packages = glob.glob(os.path.join(prefix, "lib", "python3*", "site-packages"))
    if not site_packages:
      raise RuntimeError(f"No site-packages in {prefix}")
    download(url, os.path.join(site_packages[0], relative_path))

def conda_executable():
  return os.environ.get("CONDA_EXE") or shutil.which("conda")

def conda_environments(conda):
  """Prefixes of the environments known to conda: its installation's, the user's and the registered ones."""
  try:
    result = subprocess.run([conda, "info", "--json"], check=True, capture_output=True, text=True)
    return json.loads(result.stdout).get("envs", [])
  except (OSError, subprocess.CalledProcessError, ValueError) as error:
    print(f"Could not list the conda environments ({error}), the engines will be installed.", flush=True)
    return []

def is_engine_environment(tool, prefix):
  engine = ENGINES[tool]
  if os.path.basename(os.path.normpath(prefix)) != engine["conda_environment"]:
    return False
  if not os.path.isfile(os.path.join(prefix, "bin", engine["command"])):
    return False
  if engine["distribution"]:
    name, version = engine["distribution"]
    pattern = os.path.join(prefix, "lib", "python3*", "site-packages", f"{name}-{version}.dist-info")
    return bool(glob.glob(pattern))
  return True

def find_conda_environment(tool, prefixes):
  return next((prefix for prefix in prefixes if is_engine_environment(tool, prefix)), None)

def conda_activation_env():
  env = { key: value for key, value in os.environ.items() if not key.startswith("CONDA_") or key == "CONDA_EXE" }
  env["CONDA_SHLVL"] = "0"
  return env

def install_engines(engines_dir, tools, reuse_conda=True):
  """Install the environments of the given engines and their launchers; returns {tool: launcher}.
  With reuse_conda, an existing conda environment of the engine is used instead of installing it."""
  engines_dir = os.path.abspath(os.path.expanduser(engines_dir))
  os.makedirs(os.path.join(engines_dir, "bin"), exist_ok=True)

  launchers = {}
  conda = conda_executable() if reuse_conda else None
  prefixes = conda_environments(conda) if conda else []
  for tool in tools:
    prefix = find_conda_environment(tool, prefixes)
    if prefix:
      print(f"Using the existing conda environment {prefix} for {tool}", flush=True)
      launchers[tool] = write_launcher(
        engines_dir, tool, [conda, "shell.bash", "activate", prefix], env=conda_activation_env())
  tools = [ tool for tool in tools if tool not in launchers ]
  if not tools:
    return launchers

  for name in ("pixi.toml", "pixi.lock"):
    shutil.copy2(os.path.join(manifest_dir(), name), os.path.join(engines_dir, name))
  manifest = os.path.join(engines_dir, "pixi.toml")
  pixi = pixi_binary(engines_dir)

  for tool in tools:
    environment = ENGINES[tool]["environment"]
    print(f"Installing the {tool} environment in {environment_prefix(engines_dir, tool)}", flush=True)
    run([pixi, "install", "--manifest-path", manifest, "--locked", "-e", environment], pixi_env(engines_dir))
    add_engine_files(engines_dir, tool)
    launchers[tool] = write_launcher(
      engines_dir, tool,
      [pixi, "shell-hook", "--manifest-path", manifest, "--locked", "-e", environment, "--shell", "bash"],
      env=pixi_env(engines_dir))
    for command in ENGINES[tool]["post_install"]:
      run([launchers[tool]] + command, os.environ)
  return launchers
