import os
import json
import argparse
import subprocess
import pandas as pd
from predict_binding_affinity_single_pdb import eval_binding_affinity

def load_json(json_path):
    with open(json_path) as fin:
        info = json.load(fin)
    return info

def main(input_csv, json_path, alphafold=True, save_dir=None):
    # create output folder
    if save_dir is None:
        save_dir = "./prediction_output/"
    os.makedirs(save_dir, exist_ok=True)
    # load datat
    df = pd.read_csv(input_csv)
    info = load_json(json_path)
    heavy_chain = info['heavy_chain']
    light_chain = info['light_chain']
    antigen_chains = info['antigen_chains']
    # update with predicted binding affinity value
    for idx, row in df.iterrows():
        pdb_path = row['mutated_pdb_path']
        logKd = eval_binding_affinity(pdb_path, heavy_chain, light_chain, antigen_chains, alphafold)
        df.at[idx, "predicted -log Kd"] = -logKd
    fname = os.path.splitext(os.path.basename(input_csv))[0]
    df.to_csv(f'{save_dir}/{fname}.ANTIPASTI.csv', index=False)
    return df

def parse_args():
    parser = argparse.ArgumentParser(description='Predict log(Kd) from antibody-antigen complex PDB')
    parser.add_argument('--model_output', type=str, required=True, help='Path to csv file, required column: mutated_pdb_path')
    parser.add_argument('--json_path', type=str, required=True, help='Path to PDB metadata json file')
    parser.add_argument('--no-alphafold', dest='alphafold', action='store_false',
                        help='Disable if structure is NOT from AlphaFold (default: True)')
    parser.set_defaults(alphafold=True)
    parser.add_argument('--save_dir', type=str, default="./prediction_output", help='Path to output file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.model_output, args.json_path, args.alphafold, args.save_dir)
