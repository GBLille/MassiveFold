"""
Each engine runs from its own pixi environment, (see pixi.toml and pixi.lock).
Installed when `massivefold install` is run in an engines directory with one 
launcher per engine, used as `engine_prefix` in the jobfiles (`<launcher> <command> <args>`).
"""

import glob
import hashlib
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
ALPHAFOLD3_VERSION = "1.2.0"

# scripts to download
ENGINES = {
  "AFmassive": {
    "environment": "afmassive",
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
    "scripts": {},
    "package_files": {},
    "post_install": [],
  },
  "AlphaFold3": {
    "environment": "alphafold3",
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

def pixi_env():
  env = dict(os.environ)
  env.setdefault("CONDA_OVERRIDE_CUDA", "12")
  return env

def run(command, **kwargs):
  print(" ".join(command), flush=True)
  subprocess.run(command, check=True, env=pixi_env(), **kwargs)

def write_launcher(engines_dir, tool, pixi, manifest):
  environment = ENGINES[tool]["environment"]
  activation = os.path.join(engines_dir, f"activate_{environment}.sh")
  with open(activation, "w") as handle:
    subprocess.run(
      [pixi, "shell-hook", "--manifest-path", manifest, "--frozen", "-e", environment, "--shell", "bash"],
      check=True, env=pixi_env(), stdout=handle)
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

def install_engines(engines_dir, tools):
  engines_dir = os.path.abspath(os.path.expanduser(engines_dir))
  os.makedirs(os.path.join(engines_dir, "bin"), exist_ok=True)
  for name in ("pixi.toml", "pixi.lock"):
    shutil.copy2(os.path.join(manifest_dir(), name), os.path.join(engines_dir, name))
  manifest = os.path.join(engines_dir, "pixi.toml")
  pixi = pixi_binary(engines_dir)

  launchers = {}
  for tool in tools:
    environment = ENGINES[tool]["environment"]
    print(f"Installing the {tool} environment in {environment_prefix(engines_dir, tool)}", flush=True)
    run([pixi, "install", "--manifest-path", manifest, "--frozen", "-e", environment])
    add_engine_files(engines_dir, tool)
    launchers[tool] = write_launcher(engines_dir, tool, pixi, manifest)
    for command in ENGINES[tool]["post_install"]:
      run([launchers[tool]] + command)
  return launchers
