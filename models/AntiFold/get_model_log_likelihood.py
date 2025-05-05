import os
import argparse
import pandas as pd
import numpy as np
import torch
import json
import esm
from multichain_util import extract_coords_from_complex  
from util import load_structure, get_sequence_loss  

import esm.pretrained
esm.pretrained._has_regression_weights = lambda model_name: False

def parse_args():
    parser = argparse.ArgumentParser(description="ESM-IF Benchmarking Script")
    parser.add_argument("--name", required=True, 
                        help="Key name from the JSON file (e.g., 3gbn, 4fqi, etc.)")
    parser.add_argument("--json_file", default="./data/metadata.json",
                        help="Path to the JSON file with dataset metadata")
    parser.add_argument("--gpu", type=int, default=0,
                        help="Which GPU to use; set -1 for CPU")
    parser.add_argument("--output_dir", default="./notebooks/scoring_outputs",
                        help="Directory to store output CSVs")
    parser.add_argument("--checkpoint", default="./models/AntiFold/model.pt",
                        help="Path to the model checkpoint file")
    args = parser.parse_args()
    
    return args

def score_sequence_in_complex_(model, alphabet, mutated_seqs, coords, order=None, target_chain_ids=None):
    """
    Computes the log-likelihood of the mutated complex using the ESM-IF model.
    """
    def _concatenate_coords(coords, order, padding_length=10):
        pad_coords = np.full((padding_length, 3, 3), np.nan, dtype=np.float32)
        coords_list, coords_chains = [], []
        for idx, chain_id in enumerate(order):
            if idx > 0:
                coords_list.append(pad_coords)
                coords_chains.extend(['pad'] * padding_length)
            coords_list.append(list(coords[chain_id]))
            coords_chains.extend([chain_id] * coords[chain_id].shape[0])
        coords_concatenated = np.concatenate(coords_list, axis=0)
        coords_chains = np.array(coords_chains)
        return coords_concatenated, coords_chains

    def _concatenate_seqs(mut_seqs, order, padding_length=10):
        mut_seqs_list = []
        for idx, chain_id in enumerate(order):
            if idx > 0:
                mut_seqs_list.append(['<mask>'] * (padding_length - 1) + ['<cath>'])
            mut_seqs_list.append(list(mut_seqs[chain_id]))
        mut_seqs_concatenated = ''.join([''.join(seq) for seq in mut_seqs_list])
        return mut_seqs_concatenated

    # Concatenate coordinates and sequences
    all_coords, coords_chains = _concatenate_coords(coords, order=order)
    all_seqs = _concatenate_seqs(mutated_seqs, order=order)
    # (Optional) Print debug info:
    # print("Concatenated sequence:", all_seqs)
    # print("All coordinates shape:", all_coords.shape)
    # print("Coordinates chains shape:", coords_chains.shape)
    # print("Chain order:", order)

    # Compute sequence loss using the model
    loss, target_padding_mask = get_sequence_loss(model, alphabet, all_coords, all_seqs)
    # print("Loss shape:", loss.shape)

    # Ensure dimensions match
    assert all_coords.shape[0] == coords_chains.shape[0] == loss.shape[0], "Mismatch in dimensions"

    # Compute mean loss excluding padding
    mean_loss_all = np.mean(loss[coords_chains != 'pad'])
    ll_complex = -mean_loss_all  # Negative mean loss for the complex  
    return ll_complex

def main():
    args = parse_args()

    # Load JSON metadata
    with open(args.json_file, "r") as f:
        dataset_info = json.load(f)

    # Check that the requested name is in the JSON
    if args.name not in dataset_info:
        raise ValueError(f"'{args.name}' not found in {args.json_file}.")

    metadata = dataset_info[args.name]
    
    pdb_file = metadata["pdb_path"]
    chain_order = metadata["chain_order"]
    heavy_chain = metadata["heavy_chain"]
    light_chain = metadata["light_chain"]
    antigen_chains = metadata["antigen_chains"]
    affinity_data_files = metadata["affinity_data"]
    #pdb_file = '.' + pdb_file
    #affinity_data_files = ['.' + file for file in affinity_data_files]

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)



    # 1) Create a fresh ESM-IF1 model
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()

    # 2) Load your checkpoint dictionary with torch
    ckpt = torch.load(args.checkpoint, map_location="cpu")

    # 3) Load the weights
    model.load_state_dict(ckpt)  # The entire dictionary is your param->tensor mapping
    model.eval()

    

    # If GPU is available and user requested a valid GPU ID
    if torch.cuda.is_available() and args.gpu >= 0:
        torch.cuda.set_device(args.gpu)
        model = model.cuda()
        print(f"Transferred model to GPU:{args.gpu}")
    else:
        print("Running on CPU.")

    # Load structure once (assuming the same PDB for all CSVs)
    if not os.path.exists(pdb_file):
        print(f"PDB file '{pdb_file}' does not exist. Exiting.")
        return

    try:
        structure = load_structure(pdb_file)
        coords, native_seqs = extract_coords_from_complex(structure)
        print(f"Successfully loaded structure from {pdb_file}")
        print(f"Available chain IDs in native_seqs: {list(native_seqs.keys())}")
    except Exception as e:
        print(f"Error loading/extracting from PDB '{pdb_file}': {e}")
        return

    # Loop over each CSV file specified in "affinity_data"
    for csv_path in affinity_data_files:
        if not os.path.exists(csv_path):
            print(f"CSV file '{csv_path}' not found. Skipping.")
            continue

        df = pd.read_csv(csv_path)
        print(f"Loaded {len(df)} rows from {csv_path}")

        # Prepare output file name
        csv_basename = os.path.splitext(os.path.basename(csv_path))[0]
        output_file = os.path.join(args.output_dir, f"{csv_basename}_Antifold_scores.csv")

        problematic_pdbs = []
        # Iterate over rows (variants)
        for idx, row in df.iterrows():
            # If heavy chain not in the PDB structure
            if heavy_chain not in native_seqs:
                print(f"Row {idx}: Heavy chain '{heavy_chain}' not found in PDB. Skipping.")
                problematic_pdbs.append(pdb_file)
                continue

            mutated_heavy_chain_seq = row.get('mut_heavy_chain_seq', None)
            native_heavy_seq = native_seqs[heavy_chain]
            if not mutated_heavy_chain_seq:
                print(f"Row {idx}: No 'mut_heavy_chain_seq' found in CSV. Skipping.")
                problematic_pdbs.append(pdb_file)
                continue

            if len(mutated_heavy_chain_seq) != len(native_heavy_seq):
                print(f"Row {idx}: Sequence length mismatch for heavy chain. "
                      f"Native: {len(native_heavy_seq)}, Mutated: {len(mutated_heavy_chain_seq)}")
                problematic_pdbs.append(pdb_file)
                continue

            # Build mutated sequences dict
            mutated_seqs = {heavy_chain: mutated_heavy_chain_seq}
            # Keep the native sequence for the other chains
            if light_chain in native_seqs:
                mutated_seqs[light_chain] = native_seqs[light_chain]
            else:
                print(f"Row {idx}: Light chain '{light_chain}' not found; skipping adding it.")
            for a_chain in antigen_chains:
                if a_chain in native_seqs:
                    mutated_seqs[a_chain] = native_seqs[a_chain]
                else:
                    print(f"Row {idx}: Antigen chain '{a_chain}' not found; skipping adding it.")

            # Score
            try:
                ll_complex = score_sequence_in_complex_(
                    model, alphabet,
                    mutated_seqs,
                    coords,
                    order=chain_order,
                    target_chain_ids=[heavy_chain]  # Only heavy chain is mutated
                )
            except AssertionError as e:
                print(f"Row {idx}: AssertionError during scoring: {e}")
                problematic_pdbs.append(pdb_file)
                continue
            except Exception as e:
                print(f"Row {idx}: Unexpected error during scoring: {e}")
                problematic_pdbs.append(pdb_file)
                continue

            df.at[idx, 'log-likelihood'] = ll_complex

        # Optionally drop rows missing required columns (like your original script)
        required_columns = ['log-likelihood', 'binding_score', 'mut_heavy_chain_seq']
        filtered_df = df.dropna(subset=required_columns).copy()

        # Save the results
        filtered_df.to_csv(output_file, index=False)
        print(f"Saved results to {output_file}")

        if problematic_pdbs:
            print(f"Encountered problems with: {set(problematic_pdbs)}")
        else:
            print("No problematic rows or PDB issues encountered.")

if __name__ == "__main__":
    main()