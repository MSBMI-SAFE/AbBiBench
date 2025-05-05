#!/bin/bash

# Check the number of input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 [model]"
    exit 1
fi

model=$1

# Define valid model names
valid_models=("diffab" "esmif" "dymean" "foldx" "mean" "proteinmpnn" "sasa")

# Validate the model name
if [[ ! " ${valid_models[@]} " =~ " ${model} " ]]; then
    echo "Error: Invalid model. Choose from: ${valid_models[*]}"
    exit 1
fi

# Activate the corresponding environment
eval "$(conda shell.bash hook)"
if [ "$model" = "esmif" ]; then
    conda activate struct-evo
else
    conda activate "$model"
fi
which python

# Determine the corresponding script name
script_name="run_${model}.sh"

# Check if the script exists before running it
if [ -f "$script_name" ]; then
    bash "$script_name"
else
    echo "Error: Script '$script_name' not found."
    exit 1
fi