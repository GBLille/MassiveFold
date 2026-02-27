"""Shared implementation of `massivefold run` pipeline."""

import os
import re
import glob
import json
import math
import shutil
from massivefold.parallelization import unifier
from massivefold.parallelization import batching
from massivefold.parallelization import create_jobfile

def read_json(path):
  with open(path, "r", encoding="utf-8") as handle:
    return json.load(handle)

def detect_tool_code(parameters_file):
  data = read_json(parameters_file)
  keys = list(data.keys())
  tools = [item for item in keys if item.endswith("_run")]
  if len(tools) != 1:
    raise RuntimeError("Could not detect tool section in parameter file")
  return tools[0].replace("_run", "")

def tool_from_code(short_tool):
  tool_map = {
    "AF3": "AlphaFold3",
    "AFM": "AFmassive",
    "CF": "ColabFold",
  }
  if short_tool not in tool_map:
    raise RuntimeError("Tool parameters should be one of AlphaFold3|AFmassive|ColabFold")
  return tool_map[short_tool]

def sequence_name_from_path(sequence_file):
  base = os.path.basename(sequence_file)
  if base.endswith(".fasta"):
    return base[:-6]
  return os.path.splitext(base)[0]

def next_run_name(output_dir, sequence_name, run_name):
  target = os.path.join(output_dir, sequence_name, run_name)
  if not os.path.isdir(target):
    return run_name

  print(f"Run {run_name} for {sequence_name} already exists at {target}.")
  print("Starting new iteration of this run.")
  index = 1
  candidate = run_name
  while os.path.isdir(os.path.join(output_dir, sequence_name, candidate)):
    index += 1
    candidate = f"{run_name}_{index}"
    print(f"Trying {candidate}")
  print(f"Current run is {candidate}.\n")
  return candidate

def extract_predict_times(log_dir):
  predict_times = []
  for log_file in glob.glob(os.path.join(log_dir, "jobarray_*")):
    if not os.path.isfile(log_file):
      continue
    with open(log_file, "r", encoding="utf-8") as handle:
      for line in handle:
        if "predict time" not in line:
          continue
        match = re.search(r"([0-9]*\.?[0-9]+)s\b", line)
        if match:
          predict_times.append(float(match.group(1)))
  return predict_times

def calibrated_batch_size_from_run(log_dir, wall_time, add_excess=0.1):
  predict_times = extract_predict_times(log_dir)
  if not predict_times:
    return None
  single_pred_time = math.ceil(max(predict_times))
  safe_time = math.floor(single_pred_time + (single_pred_time * add_excess))
  if safe_time <= 0:
    return None
  wall_time_seconds = wall_time * 60 * 60
  return math.floor(wall_time_seconds / safe_time)

def find_calibrated_batch_size(logs_dir, sequence_name, wall_time):
  sequence_logs_dir = os.path.join(logs_dir, sequence_name)
  if not os.path.isdir(sequence_logs_dir):
    return None

  all_runs = [
    os.path.join(sequence_logs_dir, item)
    for item in os.listdir(sequence_logs_dir)
    if os.path.isdir(os.path.join(sequence_logs_dir, item))
  ]
  if not all_runs:
    return None

  print(f"Searching for a completed preliminary run for {sequence_name}")
  lowest_pred_nb = None
  for run in all_runs:
    pred_nb = calibrated_batch_size_from_run(run, wall_time)
    if pred_nb is None:
      continue
    if lowest_pred_nb is None:
      lowest_pred_nb = pred_nb
      print(f"Found first batch size candidate: {lowest_pred_nb} at {run}")
    elif pred_nb <= lowest_pred_nb:
      lowest_pred_nb = pred_nb
      print(f"Found new batch size candidate: {lowest_pred_nb} at {run}")

  return lowest_pred_nb

def batch_size_from_specific_run(calibration_path, wall_time):
  return calibrated_batch_size_from_run(calibration_path, wall_time)

def models_from_run(path_to_run, top_n=5):
  ranking_debug = os.path.join(path_to_run, "ranking_debug.json")
  data = read_json(ranking_debug)
  models = []
  for prediction in data["order"]:
    model = prediction.split("_pred")[0]
    if model not in models:
      models.append(model)
    if len(models) == top_n:
      break
  return ",".join(models)

def copy_inputs(sequence_file, parameters_file, logs_run_dir):
  os.makedirs(logs_run_dir, exist_ok=True)
  shutil.copy2(sequence_file, logs_run_dir)
  shutil.copy2(parameters_file, logs_run_dir)

def copy_batches_file(batches_file, logs_run_dir):
  os.makedirs(logs_run_dir, exist_ok=True)
  shutil.copy2(batches_file, logs_run_dir)

def detect_precomputed_msas(tool, output_dir, sequence_name):
  if tool == "AFmassive":
    candidate = os.path.join(output_dir, sequence_name)
    if os.path.isdir(os.path.join(candidate, "msas")):
      print(f"Detected msas compatible with AFmassive for {sequence_name} at {candidate}/msas/, using them.\n")
      return candidate
  if tool == "AlphaFold3":
    candidate = os.path.join(output_dir, sequence_name, "msas_alphafold3")
    if os.path.isdir(candidate):
      print(f"Detected msas compatible with af3 for {sequence_name} at {candidate}/, using them.\n")
      return candidate
  if tool == "ColabFold":
    candidate = os.path.join(output_dir, sequence_name)
    if os.path.isdir(os.path.join(candidate, "msas_colabfold")):
      print(f"Detected msas compatible with ColabFold for {sequence_name} at {candidate}/msas_colabfold/, using them.\n")
      return candidate
  return None

def alignment_is_needed(tool, output_dir, sequence_name, force_msas_computation, msas_precomputed, waiting_for_alignment):
  if force_msas_computation:
    return True
  if waiting_for_alignment:
    return False
  if msas_precomputed:
    return False

  if tool == "AFmassive":
    return not os.path.isdir(os.path.join(output_dir, sequence_name, "msas"))
  if tool == "AlphaFold3":
    return not os.path.isdir(os.path.join(output_dir, sequence_name, "msas_alphafold3"))
  if tool == "ColabFold":
    return not os.path.isdir(os.path.join(output_dir, sequence_name, "msas_colabfold"))
  return False

def has_valid_msas(tool, msas_precomputed):
  if not msas_precomputed:
    return False
  if tool == "AFmassive":
    return os.path.isdir(os.path.join(msas_precomputed, "msas"))
  if tool == "ColabFold":
    return os.path.isdir(os.path.join(msas_precomputed, "msas_colabfold"))
  if tool == "AlphaFold3":
    return os.path.isfile(os.path.join(msas_precomputed, "msas_alphafold3_data.json"))
  return False

def safe_symlink(source, destination):
  if os.path.islink(destination):
    return
  if os.path.exists(destination):
    return
  os.symlink(source, destination)

def count_batches(batches_file):
  with open(batches_file, "r", encoding="utf-8") as handle:
    data = json.load(handle)
  return len(data)

def bool_arg(value):
  return "true" if value else "false"

def create_batches_file(
  parameters_file,
  sequence_name,
  run_name,
  predictions_per_model,
  batch_size,
  models_to_use,
  tool,
  to_screen=""):

  all_params = read_json(parameters_file)
  model_preset = batching.detect_model_preset(
    os.path.join(all_params["massivefold"]["input_dir"], f"{sequence_name}.fasta")
  )

  if tool == "AlphaFold3":
    model_names = ["AlphaFold3"]
  else:
    params_models = all_params["massivefold"].get("models_to_use", "")
    requested_models = [model.strip() for model in models_to_use.split(",") if model.strip()]
    selected_from_params = [model.strip() for model in params_models.split(",") if model.strip()]

    if model_preset == "multimer":
      model_names = [
        "model_1_multimer_v1", "model_2_multimer_v1", "model_3_multimer_v1", "model_4_multimer_v1", "model_5_multimer_v1",
        "model_1_multimer_v2", "model_2_multimer_v2", "model_3_multimer_v2", "model_4_multimer_v2", "model_5_multimer_v2",
        "model_1_multimer_v3", "model_2_multimer_v3", "model_3_multimer_v3", "model_4_multimer_v3", "model_5_multimer_v3",
      ]
    else:
      model_names = ["model_1_ptm", "model_2_ptm", "model_3_ptm", "model_4_ptm", "model_5_ptm"]

    non_existing = [model for model in selected_from_params if model not in model_names]
    non_existing.extend([model for model in requested_models if model not in model_names])
    if non_existing:
      raise RuntimeError(f"Model '{', '.join(non_existing)}' does not exist for preset '{model_preset}'")

    if selected_from_params:
      model_names = [model for model in model_names if model in selected_from_params]
    if requested_models:
      model_names = [model for model in model_names if model in requested_models]

  print(f"Running inference on models: {(', ').join(model_names)}")
  print(f"Running {predictions_per_model} predictions on each of the {len(model_names)} models")
  print(f"Total prediction number: {predictions_per_model * len(model_names)}")

  if not to_screen:
    per_model_batches = batching.batches_per_model(predictions_per_model, batch_size)
    all_model_batches = batching.batches_all_models(per_model_batches, model_names)
  else:
    all_model_batches = batching.batches_per_ligand(to_screen, predictions_per_model)

  batches_file = f"{sequence_name}_{run_name}_batches.json"
  with open(batches_file, "w", encoding="utf-8") as json_output:
    json.dump(all_model_batches, json_output, indent=4)
  return batches_file

def create_jobfile(
  job_type,
  sequence_name,
  run_name,
  path_to_parameters,
  tool,
  mf_following_msas="true",
  mf_before_inference="false",
  create_files=True):

  with open(f"{sequence_name}_{run_name}_batches.json", "r", encoding="utf-8") as json_batches:
    batches = json.load(json_batches)
  run_params = {
    "run_name": run_name,
    "sequence_name": sequence_name,
    "substitute_batch_number": list(batches.keys())[-1],
  }

  all_params = read_json(path_to_parameters)
  all_params["massivefold"]["sequence_name"] = sequence_name
  all_params["massivefold"]["run_name"] = run_name

  tool_code = "AFM" if tool == "AFmassive" else "AF3" if tool == "AlphaFold3" else "CF"
  preset_dict = {
    "model_preset": create_jobfile.detect_model_preset(
      os.path.join(all_params["massivefold"]["input_dir"], f"{sequence_name}.fasta")
    )
  }
  all_params[f"{tool_code}_run"] = preset_dict | all_params[f"{tool_code}_run"]

  run_params.update(all_params["massivefold"])
  run_params.update(all_params["custom_params"])
  run_params.update(all_params[f"{tool_code}_run"])
  run_params.update(all_params["plots"])

  if job_type == "jobarray":
    print("Parameters of the run:")
    for key in all_params["custom_params"]:
      print(f"{key}: {all_params['custom_params'][key]}")
    for key in all_params[f"{tool_code}_run"]:
      print(f"{key}: {all_params[f'{tool_code}_run'][key]}")
    print()
  if job_type == "post_treatment":
    print("Parameters for plots:")
    for key in all_params["plots"]:
      print(f"{key}: {all_params['plots'][key]}")

  if job_type != "all":
    templates = create_jobfile.group_templates(all_params, [job_type], tool)
    return create_jobfile.create_single_jobfile(
      job_type,
      templates,
      run_params,
      path_to_parameters,
      mf_following_msas,
      mf_before_inference,
      create_files,
      sequence_name,
      run_name,
    )

  templates = create_jobfile.group_templates(all_params, ["alignment", "jobarray", "post_treatment"], tool)
  return create_jobfile.create_all_jobfile(
    templates,
    run_params,
    path_to_parameters,
    mf_following_msas,
    mf_before_inference,
    create_files,
    sequence_name,
    run_name,
  )

def convert_input_if_needed(sequence_file, parameters_file, tool):
  if tool not in ["ColabFold", "AlphaFold3"]:
    return
  unifier.convert_input(sequence_file, tool, parameters_file)

def prepare_af3_inference_input(msas_json, parameters_file, batches_file):
  unifier.prepare_inference(msas_json, parameters_file, batches_file, "AlphaFold3")

def move_generated_files_to_logs(sequence_name, run_name, logs_run_dir):
  pattern = f"{sequence_name}_{run_name}_*"
  for path in glob.glob(pattern):
    target = os.path.join(logs_run_dir, os.path.basename(path))
    if os.path.exists(target):
      if os.path.isdir(target):
        shutil.rmtree(target)
      else:
        os.remove(target)
    shutil.move(path, logs_run_dir)

def submit_scheduler_job(scheduler, jobfile_content, dependency_id=None, array_size=None):
  return scheduler["submit_job"](jobfile_content, dependency_id=dependency_id, array_size=array_size)

def run_pipeline_internal(args, forwarded_args, scheduler):
  if forwarded_args:
    print("Ignoring extra arguments:", " ".join(forwarded_args))

  sequence_file = args.sequence
  run_name = args.run_name
  parameters_file = args.parameters
  predictions_per_model = args.predictions_per_model
  batch_size = args.batch_size
  wall_time = args.wall_time
  calibration = args.calibrate
  calibration_path = args.calibration_from
  msas_precomputed = args.msas_precomputed
  path_to_run = args.top_n_model
  only_msas = args.only_msas
  force_msas_computation = args.recompute_msas or args.only_msas
  wait_for_jobid = args.jobid

  if not os.path.isfile(sequence_file):
    sequence_name = sequence_name_from_path(sequence_file)
    print(f"No sequence named {sequence_name}.fasta in input directory {os.path.dirname(sequence_file)}, exiting.")
    return 1
  if not os.path.isfile(parameters_file):
    print(f"Parameter file '{parameters_file}' not found, exiting.")
    return 1

  try:
    short_tool = detect_tool_code(parameters_file)
    tool = tool_from_code(short_tool)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1

  print(f"Tool used is {tool}")
  parameters = read_json(parameters_file)
  massivefold_params = parameters.get("massivefold", {})
  output_dir = massivefold_params.get("output_dir")
  logs_dir = massivefold_params.get("logs_dir")

  if not output_dir or not logs_dir:
    print("Missing one of massivefold.output_dir|logs_dir in parameter file, exiting.")
    return 1

  sequence_name = sequence_name_from_path(sequence_file)
  run_name = next_run_name(output_dir, sequence_name, run_name)

  if not calibration and not calibration_path:
    print("No calibration for the batch size.")
  elif calibration and calibration_path:
    print("Use either -c or -C, not both, exiting.")
    return 1
  elif calibration_path and not os.path.isdir(calibration_path):
    print(f"{calibration_path} does not exist, exiting.")
    return 1
  elif calibration:
    print("Calibrating this run's batch size.")
    lowest_pred_nb = find_calibrated_batch_size(logs_dir, sequence_name, wall_time)
    if lowest_pred_nb is None:
      print(f"No preliminary run completed for {sequence_name}, exiting.")
      return 1
    batch_size = lowest_pred_nb
  elif calibration_path:
    print("Calibrating this run's batch size.")
    calibrated = batch_size_from_specific_run(calibration_path, wall_time)
    if calibrated is not None:
      batch_size = calibrated

  if calibration or calibration_path:
    print(f"Number of prediction under wall time: {batch_size}")
    if batch_size > predictions_per_model:
      batch_size = predictions_per_model
      print(f"Adapting batch size according to -p {predictions_per_model}")
    print(f"Calibrated batch size: {batch_size}\n")

  models_to_use = ""
  if path_to_run:
    print("Running with the 5 best models of the run located at path_to_run.")
    ranking_debug = os.path.join(path_to_run, "ranking_debug.json")
    if os.path.isfile(ranking_debug):
      print(f"Using {path_to_run} run to evaluate the 5 best models.")
      models_to_use = models_from_run(path_to_run)
      print(f"Using following models: {models_to_use}\n")
    else:
      print(f"Either run {path_to_run} is still running or does not exist, exiting.")
      return 1

  print(f"Run {run_name} on sequence {sequence_name} with {predictions_per_model} predictions per model")

  logs_run_dir = os.path.join(logs_dir, sequence_name, run_name)
  copy_inputs(sequence_file, parameters_file, logs_run_dir)

  batches_file = create_batches_file(
    parameters_file,
    sequence_name,
    run_name,
    predictions_per_model,
    batch_size,
    models_to_use,
    tool,
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
    force_msas_computation,
    msas_precomputed,
    waiting_for_alignment,
  ):
    convert_input_if_needed(sequence_file, parameters_file, tool)
    print(f"Running alignment for {sequence_name}")
    following_msas = not only_msas
    alignment_jobfile_content = create_jobfile(
      "alignment",
      sequence_name,
      run_name,
      parameters_file,
      tool,
      mf_following_msas=bool_arg(following_msas),
    )
    alignment_id = submit_scheduler_job(scheduler, alignment_jobfile_content)
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

    if tool == "AFmassive":
      source = os.path.realpath(os.path.join(msas_precomputed, "msas"))
      destination = os.path.join(output_dir, sequence_name, "msas")
      safe_symlink(source, destination)

    elif tool == "AlphaFold3" and not waiting_for_alignment:
      prepare_af3_inference_input(
        os.path.join(msas_precomputed, "msas_alphafold3_data.json"),
        parameters_file,
        batches_file,
      )

    shutil.copy2(parameters_file, output_run_dir)

  jobarray_jobfile_content = create_jobfile(
    "jobarray",
    sequence_name,
    run_name,
    parameters_file,
    tool,
    mf_before_inference=bool_arg(using_jobid),
  )

  array_size = count_batches(batches_file)
  dependency = alignment_id if waiting_for_alignment else None
  array_id = submit_scheduler_job(
    scheduler,
    jobarray_jobfile_content,
    dependency_id=dependency,
    array_size=array_size,
  )

  post_treatment_jobfile_content = create_jobfile("post_treatment", sequence_name, run_name, parameters_file, tool)
  submit_scheduler_job(scheduler, post_treatment_jobfile_content, dependency_id=array_id)

  move_generated_files_to_logs(sequence_name, run_name, logs_run_dir)
  return 0

def run_pipeline(args, forwarded_args, scheduler):
  try:
    return run_pipeline_internal(args, forwarded_args, scheduler)
  except RuntimeError as error:
    print(error)
    print("Exiting.")
    return 1
