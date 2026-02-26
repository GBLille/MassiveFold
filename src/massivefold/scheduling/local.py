"""Local scheduler for Linux systems."""

import itertools
import os
import subprocess
import sys

from .base import make_scheduler

_LOCAL_JOB_COUNTER = itertools.count(start=1)

def is_available():
  return sys.platform.startswith("linux")

def availability_reason():
  if is_available():
    return "Linux platform detected"
  return "Local scheduler currently supports Linux only"

def run_jobfile_once(jobfile_content, env=None):
  command = ["bash", "-c", jobfile_content]
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

def submit_job(jobfile_content, dependency_id=None, array_size=None):
  if dependency_id:
    print(f"Ignoring dependency '{dependency_id}' for local scheduler (synchronous execution)")

  if array_size is None:
    run_jobfile_once(jobfile_content)
  else:
    base_env = os.environ.copy()
    for task_id in range(array_size):
      run_jobfile_once(jobfile_content, env=array_env(base_env, task_id, array_size))

  return f"local-{next(_LOCAL_JOB_COUNTER)}"

def get_scheduler():
  return make_scheduler(
    name="local",
    description="Run locally without SLURM",
    is_available_fn=is_available,
    availability_reason_fn=availability_reason,
    submit_job_fn=submit_job,
    supports_external_dependency=False,
  )
