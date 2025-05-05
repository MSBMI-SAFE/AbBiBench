#!/usr/bin/env python3
# main.py
"""
Standard entry script: Based on the given --model and --data parameters,
first activate the corresponding conda environment (conda activate [model]),
then run the corresponding script (models/[model]/get_model_log_likelihood.py),
and pass the dataset information and other parameters to that script.
"""

import argparse
import sys
import subprocess

# Optional models
VALID_MODELS = [
    "diffab",
    "ESM-IF",
    "ESM-2",
    "dyMEAN",
    "MEAN",
    "ProteinMPNN",
    "ProSST",
    "foldx",
    "sasa"
]

# Optional datasets
VALID_DATASETS = [
    "3gbn",
    "4fqi",
    "2fjg",
    "aayl49",
    "aayl49_ml",
    "aayl51",
    "1mlc",
    "1n8z",
    "1mhp"
]

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Standard entry script for calling different model scripts on specified datasets."
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        help=f"Specify the model to use. Options are: {VALID_MODELS}"
    )
    parser.add_argument(
        "--data",
        type=str,
        required=True,
        help=f"Specify the dataset to use. Options are: {VALID_DATASETS}"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # Check if the arguments are valid
    if args.model not in VALID_MODELS:
        print(f"Error: The model '{args.model}' is not in the optional list {VALID_MODELS}.", file=sys.stderr)
        sys.exit(1)

    if args.data not in VALID_DATASETS:
        print(f"Error: The dataset '{args.data}' is not in the optional list {VALID_DATASETS}.", file=sys.stderr)
        sys.exit(1)

    # Construct the path to the script to be called
    # models/[model]/get_model_log_likelihood.py
    if args.model == 'sasa':
        script_path = f"../metrics/EpitopeSA/get_sasa.py"
        cmd = (
            f"python {script_path} --name {args.data}"
        )
    elif args.model == 'foldx':
        script_path = f"../metrics/FoldX/get_dg.py"
        cmd = (
            f"python {script_path} --name {args.data}"
        )
    else:
        script_path = f"../models/{args.model}/get_model_log_likelihood.py"
        cmd = (
            "eval \"$(conda shell.bash hook)\" "
            f"&& conda activate {args.model} "
            f"&& python {script_path} --name {args.data}"
        )

    print(cmd)


    try:
        # 'conda activate' requires shell=True to take effect
        completed_process = subprocess.run(
            cmd,
            shell=True,
            check=True,
            capture_output=False,
            text=True
        )
    except FileNotFoundError:
        print(f"Error: Could not find the script {script_path}. Please make sure the file path is correct.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print("The subprocess execution failed, return code:", e.returncode, file=sys.stderr)
        print("Subprocess standard output:", e.stdout, file=sys.stderr)
        print("Subprocess standard error:", e.stderr, file=sys.stderr)
        sys.exit(1)

    print("Subprocess script executed successfully!")
    print(f"The result is stored in: benchmark/notebooks/scoring_outputs/{args.data}_benchmarking_data_{args.model}_scores.csv")


if __name__ == "__main__":
    main()
