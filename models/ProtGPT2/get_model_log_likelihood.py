import os
import sys
import json
import argparse
import logging
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from copy import deepcopy

# Model related libraries
import transformers
from transformers import AutoTokenizer, AutoModelForCausalLM  
from Bio.SeqUtils import seq1
from torch.nn import functional as F


def parse_args():
    parser = argparse.ArgumentParser(description="ESM3 Open Full-Complex Benchmarking")
    parser.add_argument("--name", required=True,
                        help="Key name in the JSON metadata (e.g. '3gbn', '4fqi', etc.)")
    parser.add_argument("--json_file", default="./data/metadata.json",
                        help="Path to the JSON file with dataset metadata")
    parser.add_argument("--gpu", type=int, default=0,
                        help="Which GPU to use; set -1 for CPU")
    parser.add_argument("--output_dir", default="./notebooks/scoring_outputs",
                        help="Directory to store output CSVs")
    parser.add_argument("--chain_order", default='HLA',
                        help="Order of chains in the complex")
    args = parser.parse_args()
    return args

def extract_sequences_from_pdb(pdb_file: str, chains):
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

def fix_sequence(sequence):
    for i in range(len(sequence)):
        if i % 61 == 0:
            sequence = sequence[:i] + '\n' + sequence[i:]
    if sequence[-1] != '\n':
        sequence = sequence + '\n'
    full_sequence = f'<|endoftext|>{sequence}<|endoftext|>'
    
    return full_sequence

def get_ll_full_complex(model, tokenizer, sequence):
    """
    Compute the average log-likelihood over *all* residues in 'sequence'.
    No masking—just feed the entire mutated complex in and sum log-probs.
    """
    tokens = tokenizer.tokenize(fix_sequence(sequence))

    #print("this is tokens:", tokens)

    input_str = "".join(tokens)

    #print("this is input_str", input_str)
    
    inputs = tokenizer(input_str, return_tensors="pt").to(model.device)

    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits  # [batch=1, seq_len, vocab_size]
        #probs = torch.softmax(logits, dim=-1)

        # 0 -> <|endoftext|>
        # 199 -> \n (special token)
        loss = -F.nll_loss(F.log_softmax(logits[..., :-1, :].squeeze(), dim=-1), 
                          inputs['input_ids'][..., 1:].squeeze(), 
                          weight=torch.tensor([0 if i == 199 or i == 0 else 1 for i in range(tokenizer.vocab_size)],  
                                              dtype=torch.float).to(model.device))

    return loss.item()

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
        logging.info("Loading ProtGPT2 model...")

        model: AutoModelForCausalLM = AutoModelForCausalLM.from_pretrained("nferruz/ProtGPT2")
        tokenizer: AutoTokenizer = AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
        model.eval()

        if torch.cuda.is_available() and args.gpu >= 0:
            torch.cuda.set_device(args.gpu)
            model = model.cuda()
            logging.info(f"ProGPT2 model loaded onto GPU {args.gpu}")
        else:
            logging.info("Using CPU.")
    except Exception as e:
        logging.error(f"Failed loading ProGPT2 model: {e}")
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
        
        out_basename = os.path.splitext(os.path.basename(csv_path))[0]
        out_csv = os.path.join(args.output_dir, f"{out_basename}_{args.chain_order}_ProtGPT2_scores.csv")
        
        if os.path.exists(out_csv):
            results = pd.read_csv(out_csv)["log-likelihood"].tolist()
        else:
            results = [np.nan] * len(df)

        # For each variant:
        for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"Scoring {csv_path}"):
            if not np.isnan(results[idx]):
                continue
            
            mutated_heavy_chain_seq = row["mut_heavy_chain_seq"]
            # Basic check
            if not isinstance(mutated_heavy_chain_seq, str) or len(mutated_heavy_chain_seq) == 0:
                results.append(None)
                continue

            if args.chain_order == 'HLA':
                mutated_complex_seq = mutated_heavy_chain_seq + wt_lc + wt_ag
            elif args.chain_order == 'AHL':
                mutated_complex_seq = wt_ag + mutated_heavy_chain_seq + wt_lc
            elif args.chain_order == 'LAH':
                mutated_complex_seq = wt_lc + wt_ag + mutated_heavy_chain_seq
            else:
                raise ValueError(f"Invalid chain order: {args.chain_order}")

            try:
                ll_score = get_ll_full_complex(model, tokenizer, mutated_complex_seq)
            except Exception as e:
                logging.error(f"Error scoring row {idx}: {e}")
                f["log-likelihood"] = results
                df.to_csv(out_csv, index=False)
                #logging.info(f"Saved scored CSV to {out_csv}")

            results[idx] = ll_score

            if int(idx + 1) % 400 == 0:
                #print("saving")
                df["log-likelihood"] = results
                df.to_csv(out_csv, index=False)
                #logging.info(f"Saved scored CSV to {out_csv}")

        
        df["log-likelihood"] = results
        df.to_csv(out_csv, index=False)
        #logging.info(f"Saved scored CSV to {out_csv}")

    logging.info("All done!")

if __name__ == "__main__":
    main()
