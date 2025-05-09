import os
import sys
import json
import argparse
import logging
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from Bio.SeqUtils import seq1

from antiberty import AntiBERTyRunner

def parse_args():
    parser = argparse.ArgumentParser(description="AntiBERTy benchmarking")
    parser.add_argument("--name", required=True,
                        help="Key name in the JSON metadata (e.g. '3gbn', '4fqi', etc.)")
    parser.add_argument("--json_file", default="./data/metadata.json",
                        help="Path to the JSON file with dataset metadata")
    parser.add_argument("--gpu", type=int, default=0,
                        help="Which GPU to use; set -1 for CPU")
    parser.add_argument("--output_dir", default="./notebooks/scoring_outputs",
                        help="Directory to store output CSVs")
    parser.add_argument("--scoring", type=str, default="Ab-Ag", choices=["Ab", "Ab-Ag"],
                        help="Options for AntiBERTy scoring type:\n"
                             " 'Ab'    = Use PLL on heavy and light chains separately, then average. (as described in AntiBERTy GitHub)\n"
                             " 'Ab-Ag' = Use PLL on the entire complex (heavy+light+antigen).")
    args = parser.parse_args()
    return args

def extract_sequences_from_pdb(pdb_file, chains):
    """
    Extracts chain sequences from a PDB, returning (heavy_seq, light_seq, antigen_seq)
    for the chains specified in the order [heavy, light, *antigen(s)].
    """
    from Bio import PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)

    heavy_chain = chains[0] if len(chains) > 0 else None
    light_chain = chains[1] if len(chains) > 1 else None
    antigen_chains = chains[2:] if len(chains) > 2 else []

    heavy_seq, light_seq, antigen_seq = "", "", ""

    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                if residue.get_resname() in PDB.Polypeptide.standard_aa_names:
                    seq.append(seq1(residue.get_resname()))
            chain_seq = "".join(seq)

            if chain.id == heavy_chain:
                heavy_seq = chain_seq
            elif chain.id == light_chain:
                light_seq = chain_seq
            elif chain.id in antigen_chains:
                antigen_seq += chain_seq  # Concatenate multiple antigen chains

    return heavy_seq, light_seq, antigen_seq


def get_ll_antibody(antiberty, mut_hc_seq, lc_seq, batch_size = 16):
    """
    Compute pseudo log-likelihood of heavy and light chain and report their average to obtain pll of antibody as generated by pseudo_log_likelihood function from AntiBERTy.
    """

    antiberty = AntiBERTyRunner()

    sequences = [mut_hc_seq,lc_seq]

    with torch.no_grad():
        plls = antiberty.pseudo_log_likelihood(sequences, batch_size=batch_size)

    return float(plls.mean())

def get_ll_full_complex(antiberty, sequence, batch_size=16):
    """
    Compute the average log-likelihood over *all* residues in 'sequence'.
    No masking—just feed the entire mutated complex in and sum log-probs.
    """
    
    full_complex_seq = [sequence] # One combined sequence (wt_hc + lc + Ag) passed as single element list

    with torch.no_grad():
        plls = antiberty.pseudo_log_likelihood(full_complex_seq, batch_size=batch_size) 

    return float(plls[0])

def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    args = parse_args()

    if not os.path.exists(args.json_file):
        logging.error(f"JSON file '{args.json_file}' not found.")
        sys.exit(1)

    with open(args.json_file, "r") as f:
        metadata = json.load(f)

    if args.name not in metadata:
        logging.error(f"Key '{args.name}' not found in {args.json_file}")
        sys.exit(1)

    meta_info = metadata[args.name]

    pdb_path = meta_info["pdb_path"]
    heavy_chain = meta_info["heavy_chain"]
    light_chain = meta_info["light_chain"]
    antigen_chains = meta_info["antigen_chains"]
    affinity_data_files = meta_info["affinity_data"]
    # chain_order = meta_info.get("chain_order") 

    try:
        logging.info("Loading AntiBERTy model...")
        antiberty = AntiBERTyRunner()
        antiberty.model.eval()

        if torch.cuda.is_available() and args.gpu >= 0:
            torch.cuda.set_device(args.gpu)
            antiberty.model.cuda()
            logging.info(f"AntiBERTy model loaded onto GPU {args.gpu}")
        else:
            logging.info("Using CPU.")
    except Exception as e:
        logging.error(f"Failed loading AntiBERTy: {e}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(pdb_path):
        logging.error(f"PDB file '{pdb_path}' not found.")
        sys.exit(1)

    # Extract the wildtype heavy, light, and antigen sequences
    wt_hc, wt_lc, wt_ag = extract_sequences_from_pdb(
        pdb_file=pdb_path,
        chains=[heavy_chain, light_chain] + antigen_chains
    )
    logging.info(f"WT heavy: {len(wt_hc)} aa, light: {len(wt_lc)} aa, antigen(s): {len(wt_ag)} aa total")

    for csv_path in affinity_data_files:
        if not os.path.exists(csv_path):
            logging.warning(f"CSV '{csv_path}' not found. Skipping.")
            continue

        df = pd.read_csv(csv_path)
        logging.info(f"Loaded {len(df)} rows from '{csv_path}'")

        # Check for column that has the mutated heavy chain
        if "mut_heavy_chain_seq" not in df.columns:
            logging.warning(f"No 'mut_heavy_chain_seq' column in {csv_path}; skipping.")
            continue

        results = []

        # For each variant:
        for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"Scoring {csv_path}"):
            mutated_heavy_chain_seq = row["mut_heavy_chain_seq"]
            # Basic check
            if not isinstance(mutated_heavy_chain_seq, str) or len(mutated_heavy_chain_seq) == 0:
                results.append(None)
                continue

            mutated_complex_seq = mutated_heavy_chain_seq + wt_lc + wt_ag

            try:
                if args.scoring == 'Ab':
                    ll_score = get_ll_antibody(antiberty, mutated_heavy_chain_seq, wt_lc, batch_size=16)
                
                elif args.scoring == 'Ab-Ag':
                    ll_score = get_ll_full_complex(antiberty, mutated_complex_seq, batch_size=16)
            except Exception as e:
                logging.error(f"Error scoring row {idx}: {e}")
                ll_score = None

            results.append(ll_score)

        df["log-likelihood"] = results

        out_basename = os.path.splitext(os.path.basename(csv_path))[0]

        if args.scoring == 'Ab':
            out_csv = os.path.join(args.output_dir, f"{out_basename}_AntiBERTy_{args.scoring}_scores.csv")
        else:
            out_csv = os.path.join(args.output_dir, f"{out_basename}_AntiBERTy_scores.csv")
        df.to_csv(out_csv, index=False)
        logging.info(f"Saved scored CSV to {out_csv}")

    logging.info("All done!")

if __name__ == "__main__":
    main()

