#!/bin/bash

# Script to execute the Nextflow pipeline with all provided parameters
# Ensures that the specified Conda environment is activated before running Nextflow
# Usage: ./run_massivefoldNF.sh [parameters for the Nextflow pipeline]

PIPELINE_FILE="main.nf"
CONDA_ENV="nextflow"

# Check if the Conda environment exists
if ! conda info --envs | grep -q "^$CONDA_ENV"; then
    echo "Error: Conda environment '$CONDA_ENV' not found."
    exit 1
fi

# Activate the Conda environment
echo "Activating Conda environment: $CONDA_ENV"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $CONDA_ENV

# Check if the pipeline file exists
if [[ ! -f $PIPELINE_FILE ]]; then
    echo "Error: Pipeline file '$PIPELINE_FILE' not found."
    exit 1
fi

# Capture all provided arguments
NEXTFLOW_ARGS="$@"

# Run the Nextflow pipeline with the provided arguments
echo "Running: nextflow run $PIPELINE_FILE $NEXTFLOW_ARGS"
nextflow run $PIPELINE_FILE $NEXTFLOW_ARGS

# Check if the pipeline execution succeeded
if [[ $? -eq 0 ]]; then
    echo "Pipeline executed successfully."
else
    echo "Pipeline execution failed."
    exit 1
fi
