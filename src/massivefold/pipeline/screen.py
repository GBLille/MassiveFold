"""Implementation of `massivefold screen` pipeline."""

import os
import shutil

from .run import alignment_is_needed
from .run import bool_arg
from .run import copy_batches_file
from .run import copy_inputs
from .run import count_batches
from .run import create_batches_file
from .run import build_jobfile
from .run import convert_input_if_needed
from .run import detect_precomputed_msas
from .run import has_valid_msas
from .run import move_generated_files_to_logs
from .run import next_run_name
from .run import prepare_af3_inference_input
from .run import read_json
from .run import sequence_name_from_path
from .run import submit_scheduler_job

def copy_ligands_file(ligands_file, logs_run_dir):
  os.makedirs(logs_run_dir, exist_ok=True)
  shutil.copy2(ligands_file, logs_run_dir)

def run_screening_pipeline_internal(args, forwarded_args, scheduler):
  if forwarded_args:
    print("Ignoring extra arguments:", " ".join(forwarded_args))

  sequence_file = args.sequence
  ligands_file = args.ligands
  parameters_file = args.parameters
  predictions_per_model = args.predictions_per_model
  batch_size = 1
  only_msas = args.only_msas
  msas_precomputed = args.msas_precomputed
  wait_for_jobid = args.jobid
  tool = "AlphaFold3"

  if not os.path.isfile(sequence_file):
    sequence_name = sequence_name_from_path(sequence_file)
    print(f"No sequence named {sequence_name}.fasta in input directory {os.path.dirname(sequence_file)}, exiting.")
    return 1
  if not os.path.isfile(ligands_file):
    ligands_name = os.path.basename(ligands_file)
    print(f"No ligands file named {ligands_name} in input directory {os.path.dirname(ligands_file)}, exiting.")
    return 1
  if not os.path.isfile(parameters_file):
    print(f"Parameter file '{parameters_file}' not found, exiting.")
    return 1

  parameters = read_json(parameters_file)
  massivefold_params = parameters.get("massivefold", {})
  output_dir = massivefold_params.get("output_dir")
  logs_dir = massivefold_params.get("logs_dir")

  if not output_dir or not logs_dir:
    print("Missing one of massivefold.output_dir|logs_dir in parameter file, exiting.")
    return 1

  sequence_name = sequence_name_from_path(sequence_file)
  run_name = os.path.splitext(os.path.basename(ligands_file))[0]
  run_name = next_run_name(output_dir, sequence_name, run_name)

  print(f"Run {run_name} on sequence {sequence_name} with {predictions_per_model} predictions per model")

  logs_run_dir = os.path.join(logs_dir, sequence_name, run_name)
  copy_inputs(sequence_file, parameters_file, logs_run_dir)
  copy_ligands_file(ligands_file, logs_run_dir)

  batches_file = create_batches_file(
    parameters_file,
    sequence_name,
    run_name,
    predictions_per_model,
    batch_size,
    models_to_use="",
    tool=tool,
    to_screen=ligands_file,
  )
  copy_batches_file(batches_file, logs_run_dir)

  if msas_precomputed:
    print(f"Using precomputed msas at {msas_precomputed}")
  else:
    msas_precomputed = detect_precomputed_msas(tool, output_dir, sequence_name)

  waiting_for_alignment = False
  using_jobid = False
  alignment_id = None

  if wait_for_jobid:
    if not scheduler.get("supports_external_dependency", False):
      print("--jobid requires a scheduler that supports external dependencies (use --scheduler slurm).")
      return 1
    print(f"Waiting for alignment job {wait_for_jobid}")
    alignment_id = wait_for_jobid
    using_jobid = True
    waiting_for_alignment = True

  if alignment_is_needed(
    tool,
    output_dir,
    sequence_name,
    force_msas_computation=False,
    msas_precomputed=msas_precomputed,
    waiting_for_alignment=waiting_for_alignment,
  ):
    convert_input_if_needed(sequence_file, parameters_file, tool)
    print(f"Running alignment for {sequence_name}")
    alignment_jobfile_content = build_jobfile(
      "alignment",
      sequence_name,
      run_name,
      parameters_file,
      tool,
      mf_following_msas=bool_arg(not only_msas),
    )
    alignment_jobfile_name = f"alignment-{sequence_name}"
    alignment_id = submit_scheduler_job(scheduler, alignment_jobfile_content, alignment_jobfile_name)
    waiting_for_alignment = True

    if only_msas:
      move_generated_files_to_logs(sequence_name, run_name, logs_run_dir)
      print("Only run sequence alignment.")
      return 0

  elif waiting_for_alignment:
    print(f"Running inference with dependency on jobid {alignment_id}")

  elif not has_valid_msas(tool, msas_precomputed):
    print(f"Directory {msas_precomputed} does not exits or does not contain msas.")
    return 1

  else:
    print(f"{msas_precomputed} are valid.")
    output_run_dir = os.path.join(output_dir, sequence_name, run_name)
    os.makedirs(output_run_dir, exist_ok=True)

    if not waiting_for_alignment:
      prepare_af3_inference_input(
        os.path.join(msas_precomputed, "msas_alphafold3_data.json"),
        parameters_file,
        batches_file,
      )

    shutil.copy2(parameters_file, output_run_dir)

  jobarray_jobfile_content = build_jobfile(
    "jobarray",
    sequence_name,
    run_name,
    parameters_file,
    tool,
    mf_before_inference=bool_arg(using_jobid),
  )

  array_size = count_batches(batches_file)
  dependency = alignment_id if waiting_for_alignment else None

  inference_jobfile_name = f"inference-{sequence_name}_{run_name}"
  array_id = submit_scheduler_job(
    scheduler,
    jobarray_jobfile_content,
    jobfile_name=inference_jobfile_name,
    dependency_id=dependency,
    array_size=array_size,
  )

  post_treatment_jobfile_content = build_jobfile("post_treatment", sequence_name, run_name, parameters_file, tool)
  post_treatment_jobfile_name = f"post_treatment-{sequence_name}_{run_name}"
  submit_scheduler_job(scheduler, post_treatment_jobfile_content, dependency_id=array_id, jobfile_name=post_treatment_jobfile_name)

  move_generated_files_to_logs(sequence_name, run_name, logs_run_dir)
  return 0

def screening_pipeline(args, forwarded_args, scheduler):
  try:
    return run_screening_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
