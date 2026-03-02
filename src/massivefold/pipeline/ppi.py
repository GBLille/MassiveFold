"""Implementation of `massivefold screen` pipeline."""

import os
import shutil

from .screen import run_screening_pipeline_internal

def run_ppi_pipeline_internal(args, forwared_args, scheduler):
  pass

def ppi_pipeline(args, forwarded_args, scheduler):
  try:
    return run_ppi_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
