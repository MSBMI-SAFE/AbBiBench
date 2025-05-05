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


# ESM3 related libraries
import esm
from esm.sdk.api import ESM3InferenceClient, ESMProtein, ESMProteinTensor, LogitsConfig
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.structure.protein_complex import ProteinComplex


def parse_args():
    parser = argparse.ArgumentParser(description="ESM3 Large Full-Complex Benchmarking")
    parser.add_argument("--name", required=True,
                        help="Key name in the JSON metadata (e.g. '3gbn', '4fqi', etc.)")
    parser.add_argument("--json_file", default="./data/metadata.json",
                        help="Path to the JSON file with dataset metadata")
    parser.add_argument("--gpu", type=int, default=0,
                        help="Which GPU to use; set -1 for CPU")
    parser.add_argument("--output_dir", default="./notebooks/scoring_outputs",
                        help="Directory to store output CSVs")
    parser.add_argument("--with_structure", action='store_true', default=False,
                        help="Include structure in the scoring")
    args = parser.parse_args()
    return args

'''def extract_sequences_from_pdb(pdb_file: str, chains: list[str]):
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

    return heavy_seq, light_seq, antigen_seq'''

def create_wildtype(seq_dict: dict, order: str) -> ProteinComplex:
    heavy_chain_seq = seq_dict['heavy_chain_seq']
    light_chain_seq = seq_dict['light_chain_seq']
    antigen_seq = seq_dict['antigen_seqs']

    result_seq = []
    for i, seq_name in enumerate(order):
        if seq_name == "H":
            result_seq.append(heavy_chain_seq)
        elif seq_name == 'L':
            result_seq.append(light_chain_seq)
        elif seq_name == "A":
            result_seq.extend(antigen_seq)
        else:
            raise Exception(f"{seq_name} is not a valid seq character for the specified order given.")

    return ProteinComplex.from_chains(result_seq)

def get_wt_prompt(metadata: dict, 
                  order: str = "HLA") -> tuple[ESMProtein, tuple[int, int]]:
    
    antigen_chain_ids = metadata['antigen_chains']
    protein_sequences = {'heavy_chain_seq': ProteinChain.from_pdb(metadata['pdb_path'], chain_id=metadata['heavy_chain']),
                         'light_chain_seq': ProteinChain.from_pdb(metadata['pdb_path'], chain_id=metadata['light_chain']),
                         'antigen_seqs': [ProteinChain.from_pdb(metadata['pdb_path'], chain_id=chain_id) for chain_id in antigen_chain_ids]}
    
        
    wt_protein: ESMProtein = ESMProtein.from_protein_complex(create_wildtype(protein_sequences, order))
        
    return wt_protein, protein_sequences


'''def get_ll_full_complex(model, tokenizer, sequence):
    """
    Compute the average log-likelihood over *all* residues in 'sequence'.
    No masking—just feed the entire mutated complex in and sum log-probs.
    """
    tokens = tokenizer.tokenize(sequence)

    print("this is tokens:", tokens)

    input_str = "".join(tokens)

    print("this is input_str", input_str)
    
    inputs = tokenizer(input_str, return_tensors="pt").to(model.device)

    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits  # [batch=1, seq_len, vocab_size]
        probs = torch.softmax(logits, dim=-1)

    vocab = tokenizer.get_vocab()
    ll_scores = []

    for i, token in enumerate(tokens):
        if token not in vocab:
            # If unknown token -> treat as -inf
            ll_scores.append(float("-inf"))
            continue

        token_idx = vocab[token]
        p_i = probs[0, i, token_idx]
        if p_i <= 0:
            ll_scores.append(float("-inf"))
        else:
            ll_scores.append(float(torch.log(p_i)))

    return np.mean(ll_scores)'''

def mean_log_likelihood(model: ESM3InferenceClient, protein: ESMProtein, with_structure: bool = False) -> float:
    log_likelihood_scores = []
    probs_list = []
    protein_cp = deepcopy(protein)
    if not with_structure:
        protein_cp.coordinates = None
    
    # buld vocabulary mapper
    vocab = {word: i for i, word in enumerate(esm.utils.constants.esm3.SEQUENCE_VOCAB)}
    # convert protein to ProteinTensor
    protein_tensor: ESMProteinTensor = model.encode(protein_cp)
    
    # get logits #TODO: Add logitc for when the LogitOutput maxes out
    seq_logits: esm.sdk.api.LogitsOutput = model.logits(protein_tensor, LogitsConfig(sequence=True)).logits.sequence
    # convert logits to probabilities
    probs = torch.softmax(seq_logits, dim=-1)
    
    # get token selections
    for i in range(1, probs.shape[0]-1):
        if protein_tensor.sequence[i] == vocab['|']:
            probs_list.append(0.0)
            continue
        #print(i,esm.utils.constants.esm3.SEQUENCE_VOCAB[protein_tensor.sequence[i].item()],probs[i,protein_tensor.sequence[i]])
        log_likelihood_scores.append(np.log(probs[i,protein_tensor.sequence[i]]).item())
        probs_list.append(probs[i,protein_tensor.sequence[i]])
    return torch.mean(torch.tensor(log_likelihood_scores)).item(), probs_list

def get_esm3_llm(model, sequence):
    model.encode(sequence)

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

    '''
    heavy_chain = meta_info["heavy_chain"]
    light_chain = meta_info["light_chain"]
    antigen_chains = meta_info["antigen_chains"]'''
    affinity_data_files = meta_info["affinity_data"]
    # chain_order = meta_info.get("chain_order") 

    try:
        logging.info("Loading ESM3 Large model...")

        model: ESM3InferenceClient = esm.sdk.client('esm3-large-2024-03', token='1Gp0Cwdc1dahiUw2R93Q94')

    except Exception as e:
        logging.error(f"Failed loading ESM3: {e}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(pdb_path):
        logging.error(f"PDB file '{pdb_path}' not found.")
        sys.exit(1)

    # Extract the wildtype heavy, light, and antigen sequences
    wt, sequences = get_wt_prompt(meta_info)
    logging.info(f"WT heavy: {len(sequences['heavy_chain_seq'].sequence)} aa, light: {len(sequences['light_chain_seq'].sequence)} aa, antigen(s): {len(wt.sequence)} aa total")

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
        out_csv = os.path.join(args.output_dir, f"{out_basename}_ESM3-Large{'-structure' if args.with_structure else ''}_scores.csv")
        
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

            # TODO: Fix this to include structure
            mutated_wt = deepcopy(wt)
            mutated_wt.sequence = mutated_wt.sequence.replace(sequences['heavy_chain_seq'].sequence, mutated_heavy_chain_seq)
            assert len(sequences['heavy_chain_seq'].sequence) == len(mutated_heavy_chain_seq)
            assert mutated_wt.sequence != wt.sequence
            #mutated_complex_seq = "|".join([mutated_heavy_chain_seq, sequences['light_chain_seq'].sequence, "|".join([seq.sequence for seq in sequences['antigen_seqs']])])
            #mutated_complex_seq: ESMProtein = ESMProtein(sequence=mutated_complex_seq)

            try:
                #ll_score = get_ll_full_complex(model, mutated_complex_seq)
                ll_score, probs_list = mean_log_likelihood(model, mutated_wt, with_structure=args.with_structure)
            except Exception as e:
                logging.error(f"Error scoring row {idx}: {e}")
                df["log-likelihood"] = results
                df.to_csv(out_csv, index=False)
                logging.info(f"Saved scored CSV to {out_csv}")
                sys.exit(1)

            results[idx] = ll_score

            if int(idx + 1) % 10 == 0:
                #print("saving")
                df["log-likelihood"] = results
                df.to_csv(out_csv, index=False)
                logging.info(f"Saved scored CSV to {out_csv}")
                
        df["log-likelihood"] = results
        df.to_csv(out_csv, index=False)
        logging.info(f"Saved scored CSV to {out_csv}")

    logging.info("All done!")

if __name__ == "__main__":
    main()
