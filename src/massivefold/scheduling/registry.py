"""Scheduler discovery and selection logic for MassiveFold CLI."""

from .local import get_scheduler as get_local_scheduler
from .slurm import get_scheduler as get_slurm_scheduler

def scheduler_factories():
  return {
    "local": get_local_scheduler,
    "slurm": get_slurm_scheduler,
  }

def list_schedulers():
  schedulers = {}
  for name, factory in scheduler_factories().items():
    scheduler = factory()
    schedulers[name] = {
      "name": scheduler["name"],
      "description": scheduler["description"],
      "available": scheduler["is_available"](),
      "reason": scheduler["availability_reason"](),
    }
  return schedulers

def auto_select_scheduler():
  # default: slurm, fall back to local if not available 
  for candidate in ["slurm", "local"]:
    scheduler = scheduler_factories()[candidate]()
    if scheduler["is_available"]():
      return scheduler
  return scheduler_factories()["local"]()

def resolve_scheduler(requested):
  if requested == "auto":
    return auto_select_scheduler(), None

  factories = scheduler_factories()
  if requested not in factories:
    choices = ", ".join(["auto"] + sorted(factories.keys()))
    return None, f"Unknown scheduler '{requested}'. Available schedulers: {choices}."

  scheduler = factories[requested]()
  if not scheduler["is_available"]():
    return None, f"Scheduler '{requested}' is not available: {scheduler['availability_reason']()}."

  return scheduler, None
