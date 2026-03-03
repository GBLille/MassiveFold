"""Local scheduler for Linux systems."""

import itertools
import os
import re
import subprocess
import sys

from .base import make_scheduler

_LOCAL_JOB_COUNTER = itertools.count(start=1)

try:
  from tqdm import tqdm
except ImportError:
  tqdm = None

def is_available():
  return sys.platform.startswith("linux")

def availability_reason():
  if is_available():
    return "Linux platform detected"
  return "Local scheduler currently supports Linux only"

def parse_sbatch_path(jobfile_content, option_name):
  pattern_inline = re.compile(rf"^\s*#SBATCH\s+--{option_name}=([^\n]+)\s*$", re.MULTILINE)
  match = pattern_inline.search(jobfile_content)
  if match:
    return match.group(1).strip()

  pattern_spaced = re.compile(rf"^\s*#SBATCH\s+--{option_name}\s+([^\n]+)\s*$", re.MULTILINE)
  match = pattern_spaced.search(jobfile_content)
  if match:
    return match.group(1).strip()
  return None

def resolve_log_path(raw_path, job_number, task_id):
  if not raw_path:
    return None
  resolved = raw_path.replace("%j", str(job_number))
  if task_id is not None:
    resolved = resolved.replace("%a", str(task_id))
  return resolved

def run_jobfile_once(jobfile_content, env=None, log_path=None):
  command = ["bash", "-c", jobfile_content]

  if log_path and os.path.dirname(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

  if log_path:
    with open(log_path, "w", encoding="utf-8") as output_handle:
      result = subprocess.run(command, env=env, stdout=output_handle, stderr=subprocess.STDOUT)
  else:
    result = subprocess.run(command, env=env)

  if result.returncode != 0:
    raise RuntimeError(f"Job content execution failed ({result.returncode})")

def array_env(base_env, task_id, array_size):
  env = dict(base_env)
  env["SLURM_ARRAY_TASK_ID"] = str(task_id)
  env["SLURM_ARRAY_TASK_MIN"] = "0"
  env["SLURM_ARRAY_TASK_MAX"] = str(max(array_size - 1, 0))
  env["SLURM_ARRAY_TASK_COUNT"] = str(array_size)
  return env

def submit_job(jobfile_content, job_name=None, dependency_id=None, array_size=None):
  if dependency_id:
    print(f"Ignoring dependency '{dependency_id}' for local scheduler (synchronous execution)")

  job_number = next(_LOCAL_JOB_COUNTER)
  output_template = parse_sbatch_path(jobfile_content, "output")
  error_template = parse_sbatch_path(jobfile_content, "error")
  log_template = output_template or error_template
  progress_label = job_name if job_name else "local-job"
  show_progress = tqdm is not None and sys.stderr.isatty()

  if array_size is None:
    if show_progress:
      with tqdm(total=1, desc=progress_label, unit="task", leave=True) as progress:
        run_jobfile_once(
          jobfile_content,
          log_path=resolve_log_path(log_template, job_number, task_id=None),
        )
        progress.update(1)
    else:
      run_jobfile_once(
        jobfile_content,
        log_path=resolve_log_path(log_template, job_number, task_id=None),
      )
  else:
    base_env = os.environ.copy()
    if show_progress:
      with tqdm(total=array_size, desc=progress_label, unit="task", leave=True) as progress:
        for task_id in range(array_size):
          run_jobfile_once(
            jobfile_content,
            env=array_env(base_env, task_id, array_size),
            log_path=resolve_log_path(log_template, job_number, task_id=task_id),
          )
          progress.update(1)
    else:
      for task_id in range(array_size):
        run_jobfile_once(
          jobfile_content,
          env=array_env(base_env, task_id, array_size),
          log_path=resolve_log_path(log_template, job_number, task_id=task_id),
        )

  return f"local-{job_number}"

def get_scheduler():
  return make_scheduler(
    name="local",
    description="Run locally without SLURM",
    is_available_fn=is_available,
    availability_reason_fn=availability_reason,
    submit_job_fn=submit_job,
    supports_external_dependency=False,
  )
