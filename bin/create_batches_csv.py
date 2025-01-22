#!/usr/bin/env python3

import math
import json
from copy import deepcopy
from numpy import cumsum
import sys
import csv

# Choose the number of predictions inferred by each neural network model
predictions_per_model = int(sys.argv[1])
# Standard size of a prediction batch, the last batch may be smaller if the total is not a multiple
batch_size = int(sys.argv[2])
# Select models to use for prediction; None means use all models by default
models_to_use = sys.argv[3]
sequence_name =sys.argv[4]  # Name of the sequence to predict
tool = sys.argv[5] # Tool to use: "AFmassive", "AlphaFold3", or "ColabFold"

# Print the values of the variables
print(f"Predictions per model: {predictions_per_model}")
print(f"Batch size: {batch_size}")
print(f"Models to use: {models_to_use}")
print(f"Sequence name: {sequence_name}")
print(f"Tool: {tool}")
def batches_per_model(pred_nb_per_model: int, sequence_name:str):
    opt_batch_nb = math.ceil(pred_nb_per_model/batch_size)
    batch_sizes = []
    for _ in range(1, opt_batch_nb+1):
        # split total by batch of the same size, the remaining is a single smaller batch
        if pred_nb_per_model - batch_size >= 0:
            pred_nb_per_model -= batch_size
            batch_sizes.append(batch_size)
        else:
            batch_sizes.append(pred_nb_per_model)

    batch_edges = list(cumsum([0] + batch_sizes))
    one_model_batches = {i+1: {'sequence_name':sequence_name, 'start': str(batch_edges[i]), 'end': str(
        batch_edges[i+1]-1)} for i in range(len(batch_edges)-1)}
    print( one_model_batches)
    return one_model_batches

def batches_all_models(batches_unit, all_models):
    batches = {}
    for i, model in enumerate(all_models):
        unadded_batch = deepcopy(batches_unit)
        for batch in unadded_batch:
            unadded_batch[batch].update({'model': model})
            batches[str(batch+i*len(batches_unit)-1)] = unadded_batch[batch]
    return batches
            
# with open(path_to_parameters, 'r') as params:
#     all_params = json.load(params)

# if tool == "AFmassive":
#     models = all_params['massivefold']['models_to_use']
#     tool_code = "AFM"
# elif tool == "AlphaFold3":
#     models = ["AlphaFold3"]
#     tool_code = "AF3"
# elif tool == "ColabFold":
#     models = all_params['massivefold']['models_to_use']
#     tool_code = "CF"
# else:
#     raise ValueError("Unsupported tool ")

# # Get model preset and determine model names
# model_preset = all_params[f"{tool_code}_run"][f"model_preset"]

model_preset="multimer"
tool_code="CF"

if model_preset == "multimer":
    model_names = [
        f"model_{i}_{suffix}"
        for suffix in ["multimer_v1", "multimer_v2", "multimer_v3"]
        for i in range(1, 6)
    ]
elif model_preset == "monomer_ptm":
    model_names = [f"model_{i}_ptm" for i in range(1, 6)]
else:
    model_names = ["AlphaFold3"]

# Filter models if specified
if models_to_use:
    model_names = [model for model in model_names if model in models_to_use]

# Print summary
print(f"Running inference on models: {', '.join(model_names)}")
print(f"Running {predictions_per_model} predictions on each of {len(model_names)} models.")
print(f"Total prediction count: {predictions_per_model * len(model_names)}")

# Divide the predictions in batches
per_model_batches = batches_per_model( predictions_per_model, sequence_name)

# Distribute the batches on all models
all_model_batches = batches_all_models(per_model_batches, model_names)
print("bip", per_model_batches, model_names)

with open( f"{sequence_name}_batches.csv", "w") as file:
    file.write("id,sequname,start,end,model\n")
    for i in all_model_batches:
        file.write(f"{i},{sequence_name},{all_model_batches[i]['start']},{all_model_batches[i]['end']},{all_model_batches[i]['model']} \n")
    
# with open( f"{sequence_name}_batches.json", "w") as json_output:
#     json.dump(all_model_batches, json_output, indent=4)