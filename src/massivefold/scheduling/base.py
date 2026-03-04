"""Shared functional scheduler helpers for MassiveFold CLI."""

def unsupported_run(scheduler_name):
  def run(args, forwarded_args):
    raise RuntimeError(f"Scheduler '{scheduler_name}' does not implement run.")

  return run

def unsupported_screening(scheduler_name):
  def screening(args, forwarded_args):
    raise RuntimeError(f"Scheduler '{scheduler_name}' does not implement screening.")

  return screening

def unsupported_submit(scheduler_name):
  def submit(jobfile_content, dependency_id=None, array_size=None):
    raise RuntimeError(f"Scheduler '{scheduler_name}' does not implement submit_job.")

  return submit

def make_scheduler(
  name,
  description,
  is_available_fn,
  availability_reason_fn,
  run_fn=None,
  screening_fn=None,
  submit_job_fn=None,
  supports_external_dependency=False,
):
  if run_fn is None:
    run_fn = unsupported_run(name)
  if screening_fn is None:
    screening_fn = unsupported_screening(name)
  if submit_job_fn is None:
    submit_job_fn = unsupported_submit(name)

  return {
    "name": name,
    "description": description,
    "is_available": is_available_fn,
    "availability_reason": availability_reason_fn,
    "run": run_fn,
    "screening": screening_fn,
    "submit_job": submit_job_fn,
    "supports_external_dependency": supports_external_dependency,
  }
