"""SLURM scheduler using simple-slurm for submissions."""

import shutil
from simple_slurm import Slurm
from .base import make_scheduler

def is_available():
  return shutil.which("sbatch") is not None

def availability_reason():
  if shutil.which("sbatch") is None:
    return "Missing executable 'sbatch' in PATH"

def slurm_keyworded_args(job_name=None, dependency_id=None, array_size=None):
  kwargs = {}
  if job_name:
    kwargs["job_name"] = job_name
  if dependency_id:
    kwargs["dependency"] = {"afterok": dependency_id}
  if array_size is not None:
    kwargs["array"] = range(array_size)
  return kwargs

def submit_job(jobfile_content, job_name, dependency_id=None, array_size=None):
  slurm = Slurm(**slurm_keyworded_args(dependency_id=dependency_id, array_size=array_size, job_name=job_name))
  job_id = slurm.sbatch(jobfile_content)
  return str(job_id)

def get_scheduler():
  return make_scheduler(
    name="slurm",
    description="Run on a SLURM cluster",
    is_available_fn=is_available,
    availability_reason_fn=availability_reason,
    submit_job_fn=submit_job,
    supports_external_dependency=True
  )
